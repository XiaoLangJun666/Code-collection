from __future__ import print_function
import argparse
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torchvision import datasets, transforms
from torch.optim.lr_scheduler import StepLR
import matplotlib.pyplot as plt
import adabound


# Adjust the model to get a higher performance
class Net(nn.Module):
    def __init__(self):
##3 layers of 3*3
        # super(Net, self).__init__()
        # self.conv1 = nn.Conv2d(1, 5, 3, 1)  #in_channels, out_channels, kernel_size, stride=1,
        # self.conv2 = nn.Conv2d(5, 10, 3, 1)
        # self.conv3 = nn.Conv2d(10, 15, 3, 1)
        # self.dropout1 = nn.Dropout(0.25)
        # self.dropout2 = nn.Dropout(0.5)
        # self.fc1 = nn.Linear(1815, 56)
        # self.fc2 = nn.Linear(56, 10)

##2 layers of 5*5
        # super(Net, self).__init__()
        # self.conv1 = nn.Conv2d(1, 8, 5, 1)  # in_channels, out_channels, kernel_size, stride=1,
        # self.conv2 = nn.Conv2d(8, 16, 5, 1)
        # self.dropout1 = nn.Dropout(0.25)
        # self.dropout2 = nn.Dropout(0.5)
        # self.fc1 = nn.Linear(1600, 48)
        # self.fc2 = nn.Linear(48, 10)


##3 layers of convolution(3*3), 2 max pooling, 2 full connection
        # super(Net, self).__init__()
        # self.conv1 = nn.Conv2d(1, 6, 3, 1)  # in_channels, out_channels, kernel_size, stride=1,
        # self.conv2 = nn.Conv2d(6, 12, 3, 1)
        # self.conv3 = nn.Conv2d(12, 18, 3, 1)
        # self.dropout1 = nn.Dropout(0.25)
        # self.dropout2 = nn.Dropout(0.5)
        # self.fc1 = nn.Linear(450, 24)
        # self.fc2 = nn.Linear(24, 10)
##2 layers of convolution(5*5), 2 max pooling, 2 full connection
        # super(Net, self).__init__()
        # self.conv1 = nn.Conv2d(1, 8, 5, 1)  # in_channels, out_channels, kernel_size, stride=1,
        # self.conv2 = nn.Conv2d(8, 16, 3, 1)
        # self.dropout1 = nn.Dropout(0.25)
        # self.dropout2 = nn.Dropout(0.5)
        # self.fc1 = nn.Linear(256, 24)
        # self.fc2 = nn.Linear(24, 10)

##3 layers of convolution(5*5), 1 maxpooling, 2 full connection
        # super(Net, self).__init__()
        # self.conv1 = nn.Conv2d(1, 6, 5, 1)  # in_channels, out_channels, kernel_size, stride=1,
        # self.conv2 = nn.Conv2d(6, 12, 5, 1)
        # self.conv3 = nn.Conv2d(12, 18, 5, 1)
        # self.dropout1 = nn.Dropout(0.25)
        # self.dropout2 = nn.Dropout(0.5)
        # self.fc1 = nn.Linear(1152, 36)
        # self.fc2 = nn.Linear(36, 10)
##1 layer of 3*3 1 layer of 5*5 1 layer of 9*9, 1 max pooling, 2 full connection.
        # super(Net, self).__init__()
        # self.conv1 = nn.Conv2d(1, 6, 3, 1)  # in_channels, out_channels, kernel_size, stride=1,
        # self.conv2 = nn.Conv2d(6, 12, 5, 1)
        # self.conv3 = nn.Conv2d(12, 18, 9, 1)
        # self.dropout1 = nn.Dropout(0.25)  # 0.25 probability iteration failure
        # self.dropout2 = nn.Dropout(0.5)  # 0.5 probability iteration failure
        # self.fc1 = nn.Linear(882, 36)  # full connection dimension
        # self.fc2 = nn.Linear(36, 10)
##7 layer of 3*3. 1 max pooling, 2 full connection.

        super(Net, self).__init__()
        self.conv1 = nn.Conv2d(1, 6, 3, 1)  # in_channels, out_channels, kernel_size, stride=1,
        self.conv2 = nn.Conv2d(6, 12, 3, 1)
        self.conv3 = nn.Conv2d(12, 18, 3, 1)
        self.conv4 = nn.Conv2d(18, 24, 3, 1)
        self.conv5 = nn.Conv2d(24, 30, 3, 1)
        self.conv6 = nn.Conv2d(30, 36, 3, 1)
        self.conv7 = nn.Conv2d(36, 42, 3, 1)
        self.dropout1 = nn.Dropout(0.25)  # 0.25 probability iteration failure
        self.dropout2 = nn.Dropout(0.5)  # 0.5 probability iteration failure
        self.fc1 = nn.Linear(2058, 64)  # dimention of full connection
        self.fc2 = nn.Linear(64, 10)

    def forward(self, x):
        x = self.conv1(x)
        x = F.relu(x)
        x = self.conv2(x)
        x = F.relu(x)
        x = self.conv3(x)
        x = F.relu(x)
        x = self.conv4(x)
        x = F.relu(x)
        x = self.conv5(x)
        x = F.relu(x)
        x = self.conv6(x)
        x = F.relu(x)
        x = self.conv7(x)
        x = F.relu(x)
        x = F.max_pool2d(x, 2)
        x = self.dropout1(x)
        x = torch.flatten(x, 1)
        x = self.fc1(x)
        x = F.relu(x)
        x = self.dropout2(x)
        x = self.fc2(x)
        return x


def train(args, model, device, train_loader, optimizer, epoch):
    model.train()
    plt.figure()
    pic = None
    for batch_idx, (data, target) in enumerate(train_loader):
        if batch_idx in (1,2,3,4,5):
            if batch_idx == 1:
                pic = data[0,0,:,:]
            else:
                pic = torch.cat((pic,data[0,0,:,:]),dim=1)
        data, target = data.to(device), target.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = F.cross_entropy(output, target)#loss function
        # Calculate gradients
        loss.backward()
        # Optimize the parameters according to the calculated gradients
        optimizer.step()
        if batch_idx % args.log_interval == 0:
            print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                epoch, batch_idx * len(data), len(train_loader.dataset),
                100. * batch_idx / len(train_loader), loss.item()))
            if args.dry_run:
                break
    plt.imshow(pic.cpu(), cmap='gray')
    plt.show()



def test(model, device, test_loader):
    model.eval()
    test_loss = 0
    correct = 0
    with torch.no_grad():
        for data, target in test_loader:
            data, target = data.to(device), target.to(device)
            output = model(data)
            test_loss += F.cross_entropy(output, target, reduction='sum').item()  # sum up batch loss
            pred = output.argmax(dim=1, keepdim=True)  # get the index of the max log-probability
            #Calculate Confidence Values
            for i in range(0,len(pred)):
                if target[i].item()!=pred[i].item():
                    dic[str(target[i].item())] +=1
                    if str(target[i].item())=="9":
                        dic1[str(pred[i].item())]+=1

            correct += pred.eq(target.view_as(pred)).sum().item()

    test_loss /= len(test_loader.dataset)

    loss_sum.append(test_loss)
    print('\nTest set: Average loss: {:.4f}, Accuracy: {}/{} ({:.0f}%)\n'.format(
        test_loss, correct, len(test_loader.dataset),
        100. * correct / len(test_loader.dataset)))




def main():
    # Training settings
    parser = argparse.ArgumentParser(description='PyTorch 96 Example')
    parser.add_argument('--batch-size', type=int, default=64, metavar='N',
                        help='input batch size for training (default: 64)')
    parser.add_argument('--test-batch-size', type=int, default=1000, metavar='N',
                        help='input batch size for testing (default: 1000)')
    parser.add_argument('--epochs', type=int, default=14, metavar='N',
                        help='number of epochs to train (default: 14)')
    parser.add_argument('--lr', type=float, default=1, metavar='LR',
                        help='learning rate (default: 1.0)')
    parser.add_argument('--gamma', type=float, default=0.7, metavar='M',
                        help='Learning rate step gamma (default: 0.7)')
    parser.add_argument('--no-cuda', action='store_true', default=False,
                        help='disables CUDA training')
    parser.add_argument('--dry-run', action='store_true', default=False,
                        help='quickly check a single pass')
    parser.add_argument('--seed', type=int, default=1, metavar='S',
                        help='random seed (default: 1)')
    parser.add_argument('--log-interval', type=int, default=10, metavar='N',
                        help='how many batches to wait before logging training status')
    parser.add_argument('--save-model', action='store_true', default=False,
                        help='For Saving the current Model')
    args = parser.parse_args()
    use_cuda = not args.no_cuda and torch.cuda.is_available()

    torch.manual_seed(args.seed)

    device = torch.device("cuda" if use_cuda else "cpu")

    # batch_size is a crucial hyper-parameter
    train_kwargs = {'batch_size': args.batch_size}
    test_kwargs = {'batch_size': args.test_batch_size}
    if use_cuda:
        # Adjust num worker and pin memory according to your computer performance
        cuda_kwargs = {'num_workers': 1,
                       'pin_memory': True,
                       'shuffle': True}
        train_kwargs.update(cuda_kwargs)
        test_kwargs.update(cuda_kwargs)

    # Normalize the input (black and white image)
    transform=transforms.Compose([
        transforms.ToTensor(),
        transforms.Normalize((0.1307,), (0.3081,))
        ])

    # Make train dataset split
    dataset1 = datasets.MNIST('./data', train=True, download=True,
                       transform=transform)
    # Make test dataset split
    dataset2 = datasets.MNIST('./data', train=False,
                       transform=transform)

    # Convert the dataset to dataloader, including train_kwargs and test_kwargs
    train_loader = torch.utils.data.DataLoader(dataset1,**train_kwargs)
    test_loader = torch.utils.data.DataLoader(dataset2, **test_kwargs)

    # Put the model on the GPU or CPU
    model = Net().to(device)

    # Create optimizer
    optimizer = optim.Adadelta(model.parameters(), lr=args.lr) #Optimizer ：adadelta
    #optimizer = adabound.AdaBound(model.parameters(), lr=args.lr)# change Optimizer to Adabound
    # Create a schedule for the optimizer
    scheduler = StepLR(optimizer, step_size=1, gamma=args.gamma)  #调整学习率




    # Begin training and testing
    for epoch in range(1, args.epochs + 1):
        train(args, model, device, train_loader, optimizer, epoch)
        test(model, device, test_loader)
        scheduler.step()

    # Save the model
    if args.save_model:
        torch.save(model.state_dict(), "mnist_cnn.pt")




loss_sum=[]
dic={"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0}#Count the number of errors
num={"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0}#Count the total number of occurrences
dic1={"0":0,"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"8":0,"9":0}#Count the distribution of the wrong classification of the number 9
if __name__ == '__main__':
    main()
    #Loss function image
    plt.figure()
    plt.plot(loss_sum, "b--", linewidth=1)
    plt.xlabel('epoch')
    plt.ylabel("Loss")
    plt.show()
    #Distribution of statistical error classifications
    plt.figure()
    x = list(dic.keys())
    y = list(dic.values())
    plt.bar(x, y)
    plt.xlabel('class')
    plt.ylabel("error")
    plt.show()
    print(num)
    print(dic)
    #Calculate Confidence Values
    for i in dic.keys():
        a = 100 - round(dic[i] / num[i] * 100, 2)
        print("class " + str(i) + " classification confidence :" + str(a) + "%")
    print(dic1)
