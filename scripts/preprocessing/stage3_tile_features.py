import pandas as pd
import torchvision
from torchvision import transforms, datasets
import glob
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from os.path import join
from PIL import Image
import logging
import sys
import anndata as ad
import argparse
from tqdm import tqdm
import os
import pdb
import torch.nn as nn
import ResNet as ResNet
from model import Autoencoder
from tqdm import tqdm


def get_args():
    parser = argparse.ArgumentParser()

    
    '''split_index basically tells us from where we start taking the partition in dataset for that pertical job
        split_cout basically stores the  number of jobs that is created
    '''
    parser.add_argument("--tiles", dest="tiles", type=str)
    parser.add_argument("--outdir", dest="outdir", type=str)
    parser.add_argument("--num_workers", dest="num_workers", type=int, default=4)
    parser.add_argument("--batch_size", dest="batch_size", type=int, default=64)
    parser.add_argument("--njobs", dest="njobs", type=int, default=100)
    parser.add_argument("--job_i", dest="job_i", type=int, default=0)
    parser.add_argument("--debug", action="store_true", dest="debug", default=False)
    parser.add_argument("--tissue", dest="tissue", type=str)

    parser.add_argument("--dimension", dest="dimension", type=int, default=256)
    parser.add_argument("--model_type", dest="model_type", type=str)
    parser.add_argument("--simclr_ckpt", dest="simclr_ckpt", type=str, default=None)
    parser.add_argument("--autoencoder_ckpt", dest="autoencoder_ckpt", type=str, default=None)
    parser.add_argument("--kimianet_ckpt", dest="kimianet_ckpt", type=str, default=None)
    parser.add_argument("--ctranspath_ckpt", dest="ctranspath_ckpt", type=str, default=None)
    parser.add_argument("--retccl_ckpt", dest="retccl_ckpt", type=str, default=None)

    args = parser.parse_args()
    
    return args


def init_logging(logfile):
    # initialize logging
    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(stream=sys.stdout, level=logging.INFO,
        format=log_format, datefmt='%m/%d %I:%M:%S %p')
    fh = logging.FileHandler(logfile)
    fh.setFormatter(logging.Formatter(log_format))
    logging.getLogger().addHandler(fh)

class fully_connected(nn.Module):

    def __init__(self, model, num_ftrs, num_classes):
        super(fully_connected, self).__init__()
        self.model = model
        self.fc_4 = nn.Linear(num_ftrs,num_classes)

    def forward(self, x):
        x = self.model(x)
        x = torch.flatten(x, 1)
        out_1 = x
        out_3 = self.fc_4(x)
        return  out_1, out_3


def get_simclr(device, weight_path):
    from pl_bolts.models.self_supervised import SimCLR   # This import is needed only when using SIMCLR

    simclr = SimCLR.load_from_checkpoint(weight_path, strict=False)
    encoder = simclr.encoder
    encoder.fc = nn.Identity()
    encoder.eval()
    encoder = encoder.to(device)
    return encoder

def get_AutoEncoder(device, checkpoint_path):
    model = Autoencoder(1024)   # 1024 is the size of the bottleneck

    checkpoint = torch.load(checkpoint_path, map_location=device)

    # Load model weights from the checkpoint
    model.load_state_dict(checkpoint)
    # Set the model to evaluation mode
    model.eval()
    model = model.to(device)
    return model

def get_kimiaNet(device, checkpoint_path):
    
    model = torchvision.models.densenet121(pretrained=True)
    for param in model.parameters():
        param.requires_grad = False
    model.features = nn.Sequential(model.features , nn.AdaptiveAvgPool2d(output_size= (1,1)))
    num_ftrs = model.classifier.in_features
    model_final = fully_connected(model.features, num_ftrs, 30)
    model = model.to(device)
    model_final = model_final.to(device)
    model_final = nn.DataParallel(model_final)
    params_to_update = []
    criterion = nn.CrossEntropyLoss()

    model_final.load_state_dict(torch.load(checkpoint_path, map_location=device))
    
    return model

def get_ctranspath(device, checkpoint_path):
    from ctran import ctranspath
    model = ctranspath()
    model.head = nn.Identity()
    td = torch.load(checkpoint_path)
    model.load_state_dict(td['model'], strict=True)
    
    model.to(device)
    model.eval()
    
    return model

def get_retccl(device, checkpoint_path):


    
    model = ResNet.resnet50(num_classes=128,mlp=False, two_branch=False, normlinear=True)

    pretext_model = torch.load(checkpoint_path, map_location=torch.device('cpu'))
    model.fc = nn.Identity()
    model.load_state_dict(pretext_model, strict=True)


    model.eval()
    model.to(device)
    
    return model

def get_plip():
    from transformers import AutoProcessor, AutoModelForZeroShotImageClassification
    processor = AutoProcessor.from_pretrained("vinid/plip")
    model = AutoModelForZeroShotImageClassification.from_pretrained("vinid/plip")
    return model, processor

def get_model(device, tissue, model_type='retccl', simclr_ckpt=None, autoencoder_ckpt=None,
              kimianet_ckpt=None, ctranspath_ckpt=None, retccl_ckpt=None):
    
    if model_type == 'retccl':
        model = get_retccl(device, retccl_ckpt)
    elif model_type == 'simclr':
        model = get_simclr(device, simclr_ckpt)
    elif model_type == 'kimiaNet':
        model = get_kimiaNet(device, kimianet_ckpt)
    elif model_type == 'ctranspath':
        model = get_ctranspath(device, ctranspath_ckpt)
    elif model_type == 'Autoencoder':
        model = get_AutoEncoder(device, autoencoder_ckpt)
    elif model_type == 'plip':
        model, processor = get_plip()
        model.to(device)
        model.eval()
        return model, processor

    model.to(device)
    model.eval()

    return model

def imagenet_normalization():
    
    normalize = transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
    
    return normalize


class TileDataset(Dataset):

    def __init__(self, df, dimension=256):
        self.img_dim = (dimension, dimension)

        self.df = df
        mean = (0.485, 0.456, 0.406)
        std = (0.229, 0.224, 0.225)

        self.transform = transforms.Compose(
            [
                transforms.Resize(dimension),
                transforms.ToTensor(),
                transforms.Normalize(mean = mean, std = std)
            ]
        )
        
    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        img_path = self.df['path'].iloc[idx]
        img = self.transform(Image.open(img_path).convert('RGB'))
        return img, idx


# The dataset for PLIP should be in range(0,1). Since gtex already have histo image in (0,1). We don't have to do anything
class TileDataset_plip(Dataset):

    def __init__(self, df, dimension=256):
        self.img_dim = (dimension, dimension)

        self.df = df
        mean = (0.485, 0.456, 0.406)
        std = (0.229, 0.224, 0.225)

        self.transform = transforms.Compose(
            [

                transforms.Resize(dimension),
                transforms.ToTensor(),
            ]
        )
        
    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        img_path = self.df['path'].iloc[idx]
        img = self.transform(Image.open(img_path).convert('RGB'))
        return img, idx


    
def main(args):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    model_type = args.model_type
    tissue = args.tissue
    required_ckpt_by_model = {
        'retccl': ('retccl_ckpt', args.retccl_ckpt),
        'simclr': ('simclr_ckpt', args.simclr_ckpt),
        'kimiaNet': ('kimianet_ckpt', args.kimianet_ckpt),
        'ctranspath': ('ctranspath_ckpt', args.ctranspath_ckpt),
        'Autoencoder': ('autoencoder_ckpt', args.autoencoder_ckpt),
    }
    if model_type in required_ckpt_by_model:
        arg_name, ckpt_path = required_ckpt_by_model[model_type]
        if not ckpt_path:
            raise ValueError(f"--{arg_name} is required when --model_type={model_type}")

    

    if args.debug:
        pdb.set_trace()

    # outfile
    runsdir = join(args.outdir, 'runs')
    os.makedirs(runsdir, exist_ok=True)
    outfile = join(runsdir, '%.4d_%.4d.h5ad' % (args.njobs, args.job_i))

    # logging 
    logsdir = join(args.outdir, 'logs')
    os.makedirs(logsdir, exist_ok=True)
    logfile = join(logsdir, '%.4d_%.4d.h5ad' % (args.njobs, args.job_i))
    init_logging(logfile)

    # read data and split in jobs
    df = pd.read_csv(args.tiles, sep='\t', index_col=0)
    Icv = np.floor(args.njobs * np.arange(len(df)) / len(df))
    I = Icv == args.job_i
    df = df.loc[I].copy()

    if model_type == 'ctranspath':
        dataset = TileDataset(df, dimension=224)
    elif model_type == 'plip':
        dataset = TileDataset_plip(df, dimension=args.dimension)
    else:
        dataset = TileDataset(df, dimension=args.dimension)


    train_loader = torch.utils.data.DataLoader(
        dataset, batch_size=args.batch_size, shuffle=False, drop_last=False, num_workers=args.num_workers)
 
    # define model
    if model_type == 'plip':
        model, processor = get_model(
            device, tissue, model_type, args.simclr_ckpt, args.autoencoder_ckpt,
            args.kimianet_ckpt, args.ctranspath_ckpt, args.retccl_ckpt
        )
    else:
        model = get_model(
            device, tissue, model_type, args.simclr_ckpt, args.autoencoder_ckpt,
            args.kimianet_ckpt, args.ctranspath_ckpt, args.retccl_ckpt
        )

    if args.debug:
        pdb.set_trace()

    if model_type == 'kimiaNet':
        X = np.zeros([len(df), 1000], dtype=np.float32)
    elif model_type == 'ctranspath':
        X = np.zeros([len(df), 768], dtype=np.float32)
    elif model_type == 'plip':
        X = np.zeros([len(df), 512], dtype=np.float32)
    elif model_type == 'Autoencoder':
        X = np.zeros([len(df), 1024], dtype=np.float32)
    else:
        X = np.zeros([len(df), 2048], dtype=np.float32)
  


    if model_type == 'simclr':
        for batch_i, (batch, idxs) in tqdm(enumerate(train_loader)):
            batch = batch.to(device)
            outs = model(batch)[0].data.cpu().numpy()
       
            X[idxs] = outs
            logging.info(f".. analyzed {batch_i} / {len(train_loader)} batches")
    elif model_type == 'plip':
        for batch_i, (batch, idxs) in tqdm(enumerate(train_loader)):
            batch = batch.to(device)
            inputs = processor(images=batch, return_tensors="pt")
            outs = model.get_image_features(**inputs).data.cpu().numpy()
            X[idxs] = outs
            logging.info(f".. analyzed {batch_i} / {len(train_loader)} batches")
    elif model_type == 'Autoencoder':
        for batch_i, (batch, idxs) in tqdm(enumerate(train_loader)):
            batch = batch.to(device)
            outs, embed = model(batch)
            embed = embed.data.cpu().numpy()
            X[idxs] = embed
            logging.info(f".. analyzed {batch_i} / {len(train_loader)} batches")
        
    else:
        for batch_i, (batch, idxs) in tqdm(enumerate(train_loader)):
            batch = batch.to(device)
            outs = model(batch).data.cpu().numpy()
            X[idxs] = outs
            logging.info(f".. analyzed {batch_i} / {len(train_loader)} batches")


    if args.debug:
        pdb.set_trace()

    # export file
    adata = ad.AnnData(X, obs=df)
    adata.write(outfile, compression="gzip")

    for key in df.keys(): df[key] = df[key].astype(str)
    

if __name__ == "__main__":
    
    args = get_args()

    main(args)
