#!/usr/bin/env python
# coding: utf-8
#====================import packages==========================
import numpy as np
import os
#os.environ["CUDA_VISIBLE_DEVICES"] = '4'
from tqdm import tqdm
from copy import deepcopy
import matplotlib.pyplot as plt
import time
import torch
import torch.nn as nn
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torch.utils.data.sampler import RandomSampler
from imgaug import augmenters as iaa
# We try different network architecture here,including Unet、Res-Unet、ResUnetPlus，Comment or uncomment some lines to adopt
# specific network
#from Leaky_Relu import UNet
#from Res50_Leaky import UNet
#from Res18_4layer import UNet
#from core.res_unet import ResUnet
#from core.res_unet_plus import ResUnetPlusPlus
#from core.resunet_6layer import ResUnet
from Res18_Leaky import UNet
import random
batch_size=1
device = torch.device('cuda')
#define DataSet class
class TumorDataset(Dataset):
  def __init__(self, image_list, mode, mask_list=None):
    self.imagelist = image_list
    self.mode = mode
    self.masklist = mask_list
  def __len__(self):
    return len(self.imagelist)

  def __getitem__(self, idx):
    image = deepcopy(self.imagelist[idx])
    var=np.var(image)**0.5
    mean=np.mean(image)
    if mean-0<=1e-5:
        pass
    else:
        image=(image-mean)/var
    if self.mode == 'train':
      mask = deepcopy(self.masklist[idx])
      label = np.where(mask.sum() == 0, 1, 0).astype(np.float32)
      image = np.expand_dims(image, axis=0)
      mask = np.expand_dims(mask, axis=0)
      return image, mask, label

    elif self.mode == 'val':
      mask = deepcopy(self.masklist[idx])
      image = np.expand_dims(image, axis=0)
      mask = np.expand_dims(mask, axis=0)
      return image, mask

    elif self.mode == 'test':
      image = np.expand_dims(image, axis=0)
      return image


#首先加载模型
model = UNet(1)
pthfile = r'D:\学校相关申请\寒假项目\winter_school_CV_proj\要交的代码\final_model\res18_unet58.pth'
model.load_state_dict(torch.load(pthfile))
model.to(device)

resultpath=r'D:\学校相关申请\寒假项目\winter_school_CV_proj\要交的代码\result/'
datapath=r"D:\学校相关申请\寒假项目\winter_school_CV_proj\要交的代码\test\imgs"
sample=[]
all_files = os.listdir(datapath)
for file_name in all_files:
    imgpath = os.path.join(datapath, file_name)
    img=np.load(imgpath)[:,:,0]
    img = np.expand_dims(img, axis=0)  # Dataset类接受形状为（number，Height，Width）的列表
    img=TumorDataset(img,'test')
    val_loader = DataLoader(img,
                            shuffle=False,
                            batch_size=1)
    for output in val_loader:
        output=output.to(device)#模型送入了CUDA那么数据也得送入CUDA
        outputs=model(output)#尽管btach
    outputs[outputs > 0.5] = 1
    outputs[outputs <= 0.5] = 0
    outputs=outputs.detach().cpu().numpy().squeeze()
    output_name=resultpath+file_name.replace('img','seg')
    np.save(output_name,outputs)
    sample.append((imgpath,output_name))
random.shuffle(sample)
with open("pair_list.txt",'w') as f:
  f.write(str(sample))

def plot_samples(x, n=10):
  i = n
  j = 2
  plt.figure(figsize=(15, 20))
  k = 1
  idx_nums = np.random.randint(len(x), size=n)
  for idx in idx_nums:
    plt.subplot(i, j, k)
    while k % 2 != 0:
      plt.imshow(np.load(x[idx][0])[:, :, 0], cmap='gray')
      plt.xlabel("Input")
      k += 1
    plt.subplot(i, j, k)
    plt.imshow(np.load(x[idx][1])[:, :], cmap='gray')
    plt.xlabel("Ground Truth")
    k += 1
  plt.tight_layout()
  plt.show()

plot_samples(sample,n=5)











