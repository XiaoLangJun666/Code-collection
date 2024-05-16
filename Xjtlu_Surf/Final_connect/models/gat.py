import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU
from torch_geometric.nn import GATConv
from torch_geometric.nn import global_max_pool as gmp

# GAT  model
class GATNet(torch.nn.Module):
    def __init__(self, num_features_xd=78, n_output=1, num_features_xt=25,
                     n_filters=32, embed_dim=128, output_dim=128, dropout=0.2):
        super(GATNet, self).__init__()

        # graph layers
        self.gcn1 = GATConv(num_features_xd, num_features_xd, heads=10, dropout=dropout)
        self.gcn2 = GATConv(num_features_xd * 10, output_dim, dropout=dropout)
        self.fc_g1 = nn.Linear(output_dim, output_dim)

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
        self.out = nn.Linear(512, n_output)

        # activation and regularization
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, data):
        # graph input feed-forward
        x, edge_index, batch = data.x, data.edge_index, data.batch

        x = F.dropout(x, p=0.2, training=self.training)
        x = F.elu(self.gcn1(x, edge_index))
        x = F.dropout(x, p=0.2, training=self.training)
        x = self.gcn2(x, edge_index)
        x = self.relu(x)
        x = gmp(x, batch)          # global max pooling
        x = self.fc_g1(x)
        x = self.relu(x)

        # protein input feed-forward:
        target1 = data.target1
        target2 = data.target2
        embedded_xt = self.embedding_xt(target2)
        conv_xt2 = self.conv_xt_2(embedded_xt)
        conv_xt2 = self.relu(conv_xt2)
        # flatten
        xt2 = conv_xt2.view(-1, 32 * 121)
        xt2 = self.fc2_xt(xt2)

        # 1d conv layers for pssm
        # embedded_xt = self.embedding_xt(target)
        target = target1.type(torch.cuda.FloatTensor)
        conv_xt1 = self.conv_xt_1(target)
        conv_xt1 = self.relu(conv_xt1)
        # flatten
        xt1 = conv_xt1.view(-1, 544)
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
