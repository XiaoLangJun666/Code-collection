import torch
import torch.nn as nn
from core.modules import ResidualConv, Upsample


class ResUnet(nn.Module):
    def __init__(self, channel=1, filters=[64, 128, 256, 512,1024,2048]):
        super(ResUnet, self).__init__()

        self.input_layer = nn.Sequential(
            nn.Conv2d(channel, filters[0], kernel_size=(7, 7), padding=(11, 11)),
            nn.BatchNorm2d(filters[0]),
            nn.LeakyReLU(),
            nn.Conv2d(filters[0], filters[0], kernel_size=3, padding=1),
        )
        self.input_skip = nn.Sequential(
            nn.Conv2d(channel, filters[0],  kernel_size=(7, 7), padding=(11, 11))
        )

        self.residual_conv_1 = ResidualConv(filters[0], filters[1], 2, 1)
        self.residual_conv_2 = ResidualConv(filters[1], filters[2], 2, 1)
        self.add_residual_conv=ResidualConv(filters[2],filters[3],2,1)
        self.add_residual_conv_2=ResidualConv(filters[3],filters[4],2,1)
        self.add_residual_conv_3=ResidualConv(filters[4],filters[4],2,1)
        self.bridge = ResidualConv(filters[4], filters[5], 2, 1)
        self.decoder6=Decoder(2048,1024+1024,1024)
        self.decoder5=Decoder(1024,1024+1024,1024)
        self.decoder4=Decoder(1024,512+512,512)
        self.decoder3=Decoder(512,256+256,256)
        self.decoder2=Decoder(256,128+128,128)
        self.decoder1=Decoder(128,64+64,64)
        self.output_layer = nn.Sequential(
            nn.Conv2d(filters[0], 32, kernel_size=(17, 17)),
            nn.BatchNorm2d(32),
            nn.LeakyReLU(),
            nn.Conv2d(32,1,1,1),
            nn.Sigmoid(),
        )


    def forward(self, x):
        # Encode
        x1 = self.input_layer(x) + self.input_skip(x)#64,256,256
        x2 = self.residual_conv_1(x1)#128,128,128
        x3 = self.residual_conv_2(x2)#256,64,64
        x4=self.add_residual_conv(x3)#512,32,32
        x5=self.add_residual_conv_2(x4)#1024,16,16
        x6=self.add_residual_conv_3(x5)#1024,8,8
        # Bridge
        x7 = self.bridge(x6)#2048,4,4
        # Decode
        d6=self.decoder6(x7,x6)#1024,8,8
        d5=self.decoder5(d6,x5)#1024,16,16
        d4=self.decoder4(d5,x4)#512,32,32
        d3=self.decoder3(d4,x3)#256,64,64
        d2=self.decoder2(d3,x2)#128,128,128
        d1=self.decoder1(d2,x1)#64,256,256
        output = self.output_layer(d1)#1,240,240

        return output
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
if __name__ == '__main__':
    model= ResUnet(1)
    print(model)