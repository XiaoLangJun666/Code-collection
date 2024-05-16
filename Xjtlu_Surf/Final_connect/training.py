import numpy as np
import pandas as pd
import sys, os
from random import shuffle
import torch
import torch.nn as nn
from models.gat import GATNet
from models.gat_gcn import GAT_GCN
from models.gcn import GCNNet
from models.ginconv import GINConvNet
from utils import *
from torch_geometric.data import DataLoader
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import roc_curve, confusion_matrix, precision_score, recall_score, auc
import time



# training function at each epoch
def train(model, device, train_loader, optimizer, epoch):
    print('Training on {} samples...'.format(len(train_loader.dataset)))
    model.train()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    loss_accumulate = 0
    start_time = time.time()
    idx = -1
    for batch_idx, data in enumerate(train_loader):
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = loss_fn(output, data.y.view(-1, 1).float().to(device))
        loss.backward()
        optimizer.step()

        total_preds = torch.cat((total_preds, output.cpu()), 0)
        total_preds = torch.sigmoid(total_preds)
        total_labels = torch.cat((total_labels, data.y.view(-1, 1).cpu().float()), 0)
        y_pred = total_preds.flatten().tolist()
        y_true = total_labels.flatten().tolist()
        loss_accumulate += loss.item()

        if batch_idx % LOG_INTERVAL == 0 and batch_idx > 0:
            elapsed = time.time() - start_time
            print('| epoch {:3d} | {:2d}/{:2d} subsets | {:3d}/{:3d} batches '
                  '| loss {:4.3f} | auc {:4.3f} | ap {:4.3f}'.format(
                epoch, idx + 1, len(train_data), batch_idx, len(train_loader),
                       loss_accumulate / LOG_INTERVAL,
                roc_auc_score(y_true, y_pred),
                average_precision_score(y_true, y_pred)))

            total_preds = torch.Tensor()
            total_labels = torch.Tensor()
            loss_accumulate = 0
            start_time = time.time()


'''
def predicting(model, device, loader):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    print('Make prediction for {} samples...'.format(len(loader.dataset)))
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            output = model(data)
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            total_labels = torch.cat((total_labels, data.y.view(-1, 1).cpu()), 0)
    return total_labels.numpy().flatten(),total_preds.numpy().flatten()
'''


def evaluate(dataloader):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    loss_accumulate = 0
    count = 0
    with torch.no_grad():
        for idx, data in enumerate(dataloader):
            data = data.to(device)
            output = model(data)

            loss = loss_fn(total_preds, total_labels)
            loss_accumulate += loss
            count += 1

            total_preds = torch.cat((total_preds, output.cpu()), 0)
            total_preds = torch.sigmoid(total_preds)
            total_labels = torch.cat((total_labels, data.y.view(-1, 1).cpu().float()), 0)
            y_label = total_labels.flatten().tolist()
            y_pred = total_preds.flatten().tolist()
        loss = loss_accumulate / count

        fpr, tpr, thresholds = roc_curve(y_label, y_pred)

        precision = tpr / (tpr + fpr + 0.00001)
        f1 = 2 * precision * tpr / (tpr + precision + 0.00001)
        thred_optim = thresholds[5:][np.argmax(f1[5:])]
        print("optimal threshold: " + str(round(thred_optim, 3)))

        y_pred_s = [1 if i else 0 for i in (y_pred >= thred_optim)]
        auroc = auc(fpr, tpr)
        auprc = average_precision_score(y_label, y_pred)
        print("AUROC: " + str(round(auroc, 3)))
        print("AUPRC: " + str(round(auprc, 3)))

        cm1 = confusion_matrix(y_label, y_pred_s)
        print('Confusion Matrix: \n', cm1)
        print('Recall: ', round(recall_score(y_label, y_pred_s), 3))
        print('Precision: ', round(precision_score(y_label, y_pred_s), 3))

        total1 = sum(sum(cm1))
        accuracy1 = (cm1[0, 0] + cm1[1, 1]) / total1
        sensitivity1 = cm1[0, 0] / (cm1[0, 0] + cm1[0, 1])
        specificity1 = cm1[1, 1] / (cm1[1, 0] + cm1[1, 1])
        print('Accuracy: ', round(accuracy1, 3))
        print('Sensitivity: ', round(sensitivity1, 3))
        print('Specificity: ', round(specificity1, 3))

        return auroc, auprc, loss


datasets = [['kinase'][int(sys.argv[1])]]
modeling = [GINConvNet, GATNet, GAT_GCN, GCNNet][int(sys.argv[2])]
model_st = modeling.__name__

cuda_name = "cuda:0"
if len(sys.argv) > 3:
    cuda_name = ["cuda:0", "cuda:1", "cuda:2", "cuda:3", "cuda:4", "cuda:5", "cuda6", "cuda7"][int(sys.argv[3])]
print('cuda_name:', cuda_name)

TRAIN_BATCH_SIZE = 512
TEST_BATCH_SIZE = 512
LR = 0.0005
LOG_INTERVAL = 20
NUM_EPOCHS = 500

print('Learning rate: ', LR)
print('Epochs: ', NUM_EPOCHS)

# Main program: iterate over different datasets
for dataset in datasets:
    print('\nrunning on ', model_st + '_' + dataset)
    processed_data_file_train = 'data/processed/processed/' + dataset + '_grayscale_cold_protein_balance_train.pt'

    processed_data_file_test = 'data/processed/processed/' + dataset + '_grayscale_cold_protein_balance_test.pt'
    if ((not os.path.isfile(processed_data_file_train)) or (not os.path.isfile(processed_data_file_test))):
        print('please run create_data.py to prepare data in pytorch format!')
    else:
        train_data = TestbedDataset(root='data/processed', dataset=dataset + '_grayscale_cold_protein_balance_train')
        test_data = TestbedDataset(root='data/processed', dataset=dataset + '_grayscale_cold_protein_balance_test')

        # make data PyTorch mini-batch processing ready
        train_loader = DataLoader(train_data, batch_size=TRAIN_BATCH_SIZE, shuffle=True)
        test_loader = DataLoader(test_data, batch_size=TEST_BATCH_SIZE, shuffle=False)

        # training the model
        device = torch.device(cuda_name if torch.cuda.is_available() else "cpu")
        model = modeling().to(device)
        loss_fn = nn.BCEWithLogitsLoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=LR)
        patience = 0
        max_auc = 0
        model_file_name = 'model_' + model_st + '_' + dataset + '_grayscale_cold_protein_balance' + '.model'
        result_file_name = 'result_' + model_st + '_' + dataset + '_grayscale_cold_protein_balance' + '.csv'
        for epoch in range(1, NUM_EPOCHS + 1):
            epoch_start_time = time.time()
            train(model, device, train_loader, optimizer, epoch + 1)
            auroc, auprc, loss = evaluate(test_loader)
            print('predicting for valid data')
            if auroc > max_auc:
                # model_max = copy.deepcopy(model)
                torch.save(model.state_dict(), model_file_name)
                with open(result_file_name, 'w') as f:
                    f.write(','.join(map(str, [epoch, auroc, auprc, loss])))
                max_auc = auroc
                print('Validation at Epoch ' + str(epoch) + ' , AUROC: ' + str(round(auroc, 3)) + ' , AUPRC: ' + str(
                    round(auprc, 3)))
                best_epoch = epoch
                print('AUROC improved at Epoch ', best_epoch)
                patience = 0
            else:
                patience += 1
                print('No improvement since epoch', best_epoch)
                if patience >= 20:
                    print('No improvement in last 20 epochs, training stop!')
                    break
