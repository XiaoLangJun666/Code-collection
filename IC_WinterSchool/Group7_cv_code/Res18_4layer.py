import torch
import torch.nn as nn
import torchvision
from collections import OrderedDict
class Decoder(nn.Module):
    def __init__(self, in_channels, middle_channels, out_channels):
        super(Decoder, self).__init__()
        self.up = nn.ConvTranspose2d(in_channels, out_channels, kernel_size=2, stride=2)
        self.conv_relu = nn.Sequential(
            nn.Conv2d(middle_channels, out_channels, kernel_size=3, padding=1),
            nn.LeakyReLU(inplace=True)
        )

    def forward(self, x1, x2):
        x1 = self.up(x1)
        x1 = torch.cat((x1, x2), dim=1)
        x1 = self.conv_relu(x1)
        return x1


class UNet(nn.Module):
    def __init__(self, n_class=1):
        super().__init__()

        self.base_model = torchvision.models.resnet18(True)
        self.base_layers = list(self.base_model.children())
        self.layer0=nn.Conv2d(1,32,kernel_size=3,padding=1)
        self.layer1 = nn.Sequential(
            #nn.MaxPool2d(kernel_size=3, stride=2, padding=1, dilation=1, ceil_mode=False),
            nn.Conv2d(32, 64, kernel_size=(7, 7), stride=(2, 2), padding=(3, 3), bias=False),
            self.base_layers[1],
            self.base_layers[2])
        self.layer2 = nn.Sequential(*self.base_layers[3:5])
        self.layer3 = self.base_layers[5]
        self.layer4 = self.base_layers[6]

        #self.bottleneck = UNet._block(256, 512, name="bottleneck")

        self.decode3 = Decoder(256, 128+128, 128)
        self.decode2 = Decoder(128, 64+64,64)
        self.decode1 = Decoder(64, 64+64, 64)
        self.decode0 = Decoder(64, 32 + 32, 32)
        # self.decode0 = nn.Sequential(
        #     # nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True),
        #     # nn.Conv2d(64, 32, kernel_size=3, padding=1, bias=False),
        #     # nn.Conv2d(32, 64, kernel_size=3, padding=1, bias=False)
        #     nn.ConvTranspose2d(64,32,kernel_size=2, stride=2)
        # )
        #self.conv_last = nn.Sequential(nn.Conv2d(64, n_class, kernel_size=(17,17)),nn.Softmax(dim=1))
        self.conv_last = nn.Sequential(nn.Conv2d(32, n_class, 1), nn.Sigmoid())
        #self.conv_last = nn.Conv2d(64, n_class, kernel_size=(17, 17))
    def forward(self, input):
        e0=self.layer0(input)#32,240,240
        e1 = self.layer1(e0)  # 64,120,120
        e2 = self.layer2(e1)  # 64,60,60
        e3 = self.layer3(e2)  # 128,30,30
        f = self.layer4(e3)  # 256,15,15
        #bottle=self.bottleneck(e4)#512,15,15
        d3 = self.decode3(f, e3)  # 128,30,30
        d2 = self.decode2(d3, e2)  # 64,60,60
        d1 = self.decode1(d2, e1)  # 64,120,120
        d0 = self.decode0(d1, e0)  # 32,240,240
        #d0 = self.decode0(d1)  # 64,256,256
        out = self.conv_last(d0)  # 1,240,240
        return out


if __name__ == '__main__':
    model=UNet()
    print(model)