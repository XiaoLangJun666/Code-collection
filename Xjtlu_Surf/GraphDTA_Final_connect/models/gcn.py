import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_max_pool as gmp
import numpy as np


# GCN based model
class GCNNet(torch.nn.Module):
    def __init__(self, n_output=1, n_filters=32, embed_dim=128,num_features_xd=78, num_features_xt=25, output_dim=128, dropout=0.2):

        super(GCNNet, self).__init__()

        # SMILES graph branch
        self.n_output = n_output
        self.conv1 = GCNConv(num_features_xd, num_features_xd)
        self.conv2 = GCNConv(num_features_xd, num_features_xd*2)
        self.conv3 = GCNConv(num_features_xd*2, num_features_xd * 4)
        self.fc_g1 = torch.nn.Linear(num_features_xd*4, 1024)
        self.fc_g2 = torch.nn.Linear(1024, output_dim)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout)

        # protein sequence branch (1d conv)for one hot
        self.embedding_xt = nn.Embedding(num_features_xt + 1, embed_dim)
        self.conv_xt_2 = nn.Conv1d(in_channels=750, out_channels=n_filters, kernel_size=8)
        self.fc2_xt = nn.Linear(32 * 121, output_dim)

        # protein sequence branch (1d conv)for pssm
        # self.embedding_xt = nn.Embedding(num_features_xt + 1, embed_dim)
        self.conv_xt_1 = nn.Conv1d(in_channels=750, out_channels=n_filters, kernel_size=8)
        self.fc1_xt = nn.Linear(416, output_dim)

        # combined layers
        self.fc1 = nn.Linear(3 * output_dim, 1024)
        self.fc2 = nn.Linear(1024, 512)
        self.out = nn.Linear(512, self.n_output)


        """
        以上部分保留了源代码onehot的embedding层之外，新设了新特征输入的卷积层，两个同时保留，注意，当我们更改输入的时候，著需要改动in_channels(输入的特征的行数),
        544(计算数据输入后得到的剩余维度，总数/512)，其余网络依次按此更改，这里的self.fc1 = nn.Linear(3 * output_dim, 1024)因为输入增加所以改为3*
        """
    def forward(self, data):
        # get graph input
        x, edge_index, batch = data.x, data.edge_index, data.batch
        # get protein input
        #indices = indices.to(torch.int64)



        target1=data.target1
        target2=data.target2

        x = self.conv1(x, edge_index)
        x = self.relu(x)

        x = self.conv2(x, edge_index)
        x = self.relu(x)

        x = self.conv3(x, edge_index)
        x = self.relu(x)
        x = gmp(x, batch)       # global max pooling

        # flatten
        x = self.relu(self.fc_g1(x))
        x = self.dropout(x)
        x = self.fc_g2(x)
        x = self.dropout(x)

        # 1d conv layers for one hot
        embedded_xt = self.embedding_xt(target2)
        conv_xt2 = self.conv_xt_2(embedded_xt)
        # flatten
        xt2 = conv_xt2.view(-1, 32 * 121)
        xt2 = self.fc2_xt(xt2)

        # 1d conv layers for pssm
        # embedded_xt = self.embedding_xt(target)
        target1 = target1.type(torch.cuda.FloatTensor)
        conv_xt1 = self.conv_xt_1(target1)
        # flatten
        xt1 = conv_xt1.view(-1, 416)
        xt1 = self.fc1_xt(xt1)

        # concat

        xc = torch.cat((x, xt1,xt2), 1)
        # add some dense layers
        xc = self.fc1(xc)
        xc = self.relu(xc)
        xc = self.dropout(xc)
        xc = self.fc2(xc)
        xc = self.relu(xc)
        xc = self.dropout(xc)
        out = self.out(xc)
        return out
