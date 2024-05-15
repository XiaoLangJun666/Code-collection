import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU
from torch_geometric.nn import GCNConv, GATConv, GINConv, global_add_pool
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp

# GCN-CNN based model

class GAT_GCN(torch.nn.Module):
    def __init__(self, n_output=1, num_features_xd=78, num_features_xt=25,
                 n_filters=32, embed_dim=128, output_dim=128, dropout=0.2):

        super(GAT_GCN, self).__init__()

        self.n_output = n_output
        self.conv1 = GATConv(num_features_xd, num_features_xd, heads=10)
        self.conv2 = GCNConv(num_features_xd*10, num_features_xd*10)
        self.fc_g1 = torch.nn.Linear(num_features_xd*10*2, 1500)
        self.fc_g2 = torch.nn.Linear(1500, output_dim)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout)

        # protein sequence branch (1d conv)for one hot
        self.embedding_xt = nn.Embedding(num_features_xt + 1, embed_dim)
        self.conv_xt_2 = nn.Conv1d(in_channels=1000, out_channels=n_filters, kernel_size=8)
        self.fc2_xt = nn.Linear(32 * 121, output_dim)

        # protein sequence branch (1d conv)for pssm
        # self.embedding_xt = nn.Embedding(num_features_xt + 1, embed_dim)
        self.conv_xt_1 = nn.Conv1d(in_channels=20, out_channels=n_filters, kernel_size=8)
        self.fc1_xt = nn.Linear(544, output_dim)

        # combined layers
        self.fc1 = nn.Linear(3 * output_dim, 1024)
        self.fc2 = nn.Linear(1024, 512)
        self.out = nn.Linear(512, self.n_output)  # n_output = 1 for regression task

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        target1 = data.target1
        target2 = data.target2
        # print('x shape = ', x.shape)
        x = self.conv1(x, edge_index)
        x = self.relu(x)
        x = self.conv2(x, edge_index)
        x = self.relu(x)
        # apply global max pooling (gmp) and global mean pooling (gap)
        x = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        x = self.relu(self.fc_g1(x))
        x = self.dropout(x)
        x = self.fc_g2(x)

        # 1d conv layers for one hot
        embedded_xt = self.embedding_xt(target2)
        conv_xt2 = self.conv_xt_2(embedded_xt)
        # flatten
        xt2 = conv_xt2.view(-1, 32 * 121)
        xt2 = self.fc2_xt(xt2)

        # 1d conv layers for pssm
        # embedded_xt = self.embedding_xt(target)
        target = target1.type(torch.cuda.FloatTensor)
        conv_xt1 = self.conv_xt_1(target)
        # flatten
        xt1 = conv_xt1.view(-1, 544)
        xt1 = self.fc1_xt(xt1)

        # concat
        xc = torch.cat((x, xt1,xt2, 1))
        # add some dense layers
        xc = self.fc1(xc)
        xc = self.relu(xc)
        xc = self.dropout(xc)
        xc = self.fc2(xc)
        xc = self.relu(xc)
        xc = self.dropout(xc)
        out = self.out(xc)
        return out
