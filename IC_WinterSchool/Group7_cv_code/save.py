#!/usr/bin/env python
# coding: utf-8
#====================import packages==========================
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt

import random


batch_size =6
epoch = 4
smooth = 1e-5



#此函数用于将得到的图片路径全部收集起来存成.npy文件，避免每次都要遍历路径
def ImageFetch(img_list,seg_list):
    image_total = np.zeros((len(img_list), 240, 240,3), dtype=np.float32)
    mask_total = np.zeros((len(img_list), 240, 240), dtype=int)

    for idx, image_id in tqdm(enumerate(zip(img_list,seg_list)), total=len(img_list)):
        image_path = image_id[0]
        mask_path = image_id[1]

        image = np.load(image_path)
        mask = np.load(mask_path)

        image_total[idx] = image
        mask_total[idx] = mask
    np.save(r'F:/Yesimg_all.npy',image_total)
    np.save(r'F:/Yesmask_all.npy',mask_total)
    return image_total, mask_total

img_path=r"D:\学校相关申请\寒假项目\winter_school_CV_proj\Dataset\Yes"
img_list=[]
mask_list=[]
all_files = os.listdir(img_path )
files = [item for item in all_files if "img" in item]
random.shuffle(files)
img_num = len(files)
for (n, file_name) in enumerate(files):
    img = os.path.join(img_path,file_name)
    seg = os.path.join(img_path,file_name.split('_')[0]+'_seg.npy')
    img_list.append(img)
    mask_list.append(seg)

# for CLASS in os.listdir(img_path):
#     if not CLASS.startswith('.'):
#         all_files = os.listdir(img_path + CLASS)
#         files = [item for item in all_files if "img" in item]
#         random.shuffle(files)
#         img_num = len(files)
#         for (n, file_name) in enumerate(files):
#             img = os.path.join(img_path,CLASS,file_name)
#             seg = os.path.join(img_path,CLASS,file_name.split('_')[0]+'_seg.npy')
#             mask_list.append(seg)
#             img_list.append(img)
# total=list(zip(img_list,mask_list))
# random.shuffle(total)
# img_list=[item[0] for item in total]
# mask_list=[item[1] for item in total]

#通过取出部分图片来检查是否保存有误
ImageFetch(img_list,mask_list)
image=np.load(r'F:/Yesimg_all.npy')
mask=np.load(r'F:/Yesmask_all.npy')
print("shape_img",np.shape(image))
print("shape_mask",np.shape(mask))
k=1
count=0
for i in range(3):
    plt.subplot(3, 4, k)
    r = image[count]
    plt.imshow(r[:, :, 0], cmap='gray')
    plt.xlabel("0")
    k += 1
    plt.subplot(3, 4, k)
    plt.imshow(r[:, :, 1], cmap='gray')
    plt.xlabel("1")
    k += 1
    plt.subplot(3, 4, k)
    plt.imshow(r[:, :, 2], cmap='gray')
    plt.xlabel("2")
    k += 1
    plt.subplot(3, 4, k)
    plt.imshow(mask[count], cmap='gray')
    plt.xlabel("truth")
    k+=1
    if k == 13:
        break
    count = count + 1
plt.tight_layout()
plt.show()
