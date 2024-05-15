#!/usr/bin/env python
# coding: utf-8
#====================import packages==========================
import numpy as np
import os
os.environ["CUDA_VISIBLE_DEVICES"] = '4'
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
## Hyperparameters
batch_size =60
epoch = 4
smooth = 1e-5
device = torch.device('cuda')
save_weight = '/data/lpy_data/final_weight'# saved path of model
if not os.path.isdir(save_weight):
  os.mkdir(save_weight)
weight_name = 'model_' + 'res18_unet'


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



# ## Train your model here
# Once you defined the model and data generator, you can start training your model.

#define loss
def dice_coef(predict,target):
    epsilon = 1e-5
    assert predict.size() == target.size(), "the size of predict and target must be equal."
    num = predict.size(0)
    pre = predict.view(num, -1)
    tar = target.view(num, -1)

    intersection = (pre * tar).sum(-1)  # 利用预测值与标签相乘当作交集
    union = (pre + tar).sum(-1)
    dice_score = (2 * intersection + epsilon) / (union + epsilon)
    return dice_score.sum()

class DiceLoss(nn.Module):
    def __init__(self):
        super(DiceLoss, self).__init__()
        self.epsilon = 1e-5

    def forward(self, predict, target):
        assert predict.size() == target.size(), "the size of predict and target must be equal."
        num = predict.size(0)
        pre = predict.view(num, -1)
        tar = target.view(num, -1)
        intersection = (pre * tar).sum(-1) # 利用预测值与标签相乘当作交集
        union = (pre + tar).sum(-1)
        dice_score=(2 * intersection + self.epsilon) / (union + self.epsilon)
        return (1 - dice_score.sum()/num)

# Load data
#X_train, y_train = ImageFetch(img_list,seg_list,'train')
#I divide the given dataset inta train and val set.And save them as .npy
X_train=np.load("/data/lpy_data/aug_train_image_total.npy")[:,:,:,0]#   F:/train_image_total.npy
y_train=np.load("/data/lpy_data/aug_train_mask_total.npy")#    F:/train_mask_total.npy
train_data = TumorDataset(X_train, 'train', y_train)
train_loader = DataLoader(train_data,
                          shuffle=RandomSampler(train_data),
                          batch_size=batch_size)

def train(train_loader, model,optimizer,loss_class):
  running_loss = 0.0
  data_size = np.shape(y_train)[0]
  model.train()
  '''
  for inputs_, masks_, labels in tqdm(train_loader):
    seq = iaa.Sometimes(0.2,
        iaa.SomeOf((0, 2), [  # 建立一个名为seq的实例，定义增强方法，用于增强
          iaa.Fliplr(1.0),
          iaa.Flipud(1.0),
          iaa.GaussianBlur(sigma=(0, 3.0)),  # 这个只会对原图操作，不会对mask操作
          iaa.Affine(translate_percent=(0, 0.1), rotate=(-15, 15), cval=0, mode='constant'),
          iaa.CropAndPad(px=(-20, 0), percent=None, pad_mode='constant', pad_cval=0, keep_size=True)
    ], random_order=True))
    images_ = np.expand_dims(inputs_.squeeze(), axis=3)
    segs_ = np.expand_dims(masks_.squeeze(), axis=3)
    images_aug, masks_aug = seq(images=images_, segmentation_maps=segs_)
    images_aug, masks_aug= torch.tensor(images_aug), torch.tensor(masks_aug)
    inputs = torch.unsqueeze(images_aug.squeeze(), dim=1).to(device)
    masks = torch.unsqueeze(masks_aug.squeeze(),dim=1).to(device)
    Here is the appraoch of data augmentation
  '''
  for inputs, masks, labels in tqdm(train_loader):
    inputs, masks, labels = inputs.to(device),masks.to(device), labels.to(device)
    with torch.set_grad_enabled(True):
      logit = model(inputs)
      #loss = nn.BCEWithLogitsLoss()(logit.squeeze(1), masks.squeeze(1).float())
      loss=loss_class(logit.squeeze(1),masks.squeeze(1))
      loss.backward()
      optimizer.zero_grad()
      optimizer.step()
    running_loss += loss.item() * inputs.size(0)
  epoch_loss = running_loss / data_size
  return epoch_loss

# ## Run the model on the test set
# After your last Q&A session, you will be given the test set. Run your model on the test set to get the segmentation results and submit your results in a .zip file. If the MRI image is named '100_img.npy', save your segmentation result as '100_seg.npy'.

# In[ ]:

X_val=np.load("/data/lpy_data/test_image_total.npy")[:,:,:,0]
y_val=np.load("/data/lpy_data/test_mask_total.npy")
val_data = TumorDataset(X_val, 'val', y_val)
val_loader = DataLoader(val_data,
                        shuffle=False,
                        batch_size=batch_size)
def test(test_loader, model,loss_class):
  running_loss = 0.0
  data_size = np.shape(y_val)[0]
  dice_total=0.0
  model.eval()
  for inputs, masks in tqdm(test_loader):
    inputs= inputs.to(device)
    masks=masks.to(device)
    with torch.no_grad():
      outputs= model(inputs)
      loss=loss_class(outputs.squeeze(1),masks.squeeze(1))
      outputs[outputs>0.5]=1
      outputs[outputs<=0.5]=0
      # for i in range(np.shape(outputs)[0]):
      #   super_threshold_indices = outputs[i] > outputs[i].max() * 0.75
      #   outputs[i][super_threshold_indices] = 1
      #   lower_threshold_indices = outputs[i] <= outputs[i].max() * 0.75
      #   outputs[i][lower_threshold_indices] = 0
      #if BCEwithLlogitsloss is adopted, uncomment above code
      dice=dice_coef(outputs,masks)
      dice_total+=dice
      running_loss += loss.item() * inputs.size(0)
  dice_mean=dice_total/data_size
  epoch_loss = running_loss / data_size
  return epoch_loss, dice_mean.detach().cpu()

#####################################################
#every epoch,a test will be done
model = UNet(1)
#model=ResUnetPlusPlus(1)
#model=ResUnet(1)
model.to(device)
# Setup optimizer
optimizer = torch.optim.Adam(model.parameters(),lr=0.0001)
loss=DiceLoss()
best_acc = 0
trainlosslist=[]
vallosslist=[]
dice_score=[]
for epoch_ in range(epoch):
    train_loss = train(train_loader, model,optimizer,loss)
    val_loss, dice = test(val_loader, model,loss)
    trainlosslist.append(train_loss)
    vallosslist.append(val_loss)
    dice_score.append(dice)

    if dice > best_acc:
        best_acc = dice
        best_param = model.state_dict()
        torch.save(best_param, save_weight + weight_name + str(epoch_)+ '.pth')
    print('epoch: {} train_loss: {:.12f} val_loss: {:.12f} val_dice_score: {:.12f}'.format(epoch_ + 1, train_loss, val_loss, dice))
print("train_loss",trainlosslist)
print("val_loss",vallosslist)
print("dice_score",dice_score)


plt.ylabe("dice")
plt.xlabel("epoch")
plt.plot(dice_score)
plt.savefig("./dice.png")
plt.clf()
plt.ylabe("train_loss")
plt.xlabel("epoch")
plt.plot(trainlosslist)
plt.savefig("./train_loss.png")
plt.clf()
plt.ylabe("val_loss")
plt.xlabel("epoch")
plt.plot(vallosslist)
plt.savefig("./val_loss.png")
print("best_dice_score",best_acc)






