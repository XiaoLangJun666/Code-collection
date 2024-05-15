import torch
import torch.nn as nn
from core.modules import ResidualConv, Upsample


class ResUnet(nn.Module):
    def __init__(self, channel=1, filters=[64, 128, 256, 512,1024]):
        super(ResUnet, self).__init__()

        self.input_layer = nn.Sequential(
            nn.Conv2d(channel, filters[0], kernel_size=3, padding=1),
            nn.BatchNorm2d(filters[0]),
            nn.LeakyReLU(),
            nn.Conv2d(filters[0], filters[0], kernel_size=3, padding=1),
        )
        self.input_skip = nn.Sequential(
            nn.Conv2d(channel, filters[0], kernel_size=3, padding=1)
        )

        self.residual_conv_1 = ResidualConv(filters[0], filters[1], 2, 1)
        self.residual_conv_2 = ResidualConv(filters[1], filters[2], 2, 1)
        self.add_residual_conv=ResidualConv(filters[2],filters[3],2,1)

        self.bridge = ResidualConv(filters[3], filters[4], 2, 1)

        self.add_upsample = Upsample(filters[4], filters[4], 2, 2)
        self.add_up_residual=ResidualConv(filters[4] + filters[3], filters[3], 1, 1)
        self.upsample_1=Upsample(filters[3], filters[3], 2, 2)
        self.up_residual_conv1 = ResidualConv(filters[3] + filters[2], filters[2], 1, 1)

        self.upsample_2 = Upsample(filters[2], filters[2], 2, 2)
        self.up_residual_conv2 = ResidualConv(filters[2] + filters[1], filters[1], 1, 1)

        self.upsample_3 = Upsample(filters[1], filters[1], 2, 2)
        self.up_residual_conv3 = ResidualConv(filters[1] + filters[0], filters[0], 1, 1)

        self.output_layer = nn.Sequential(
            nn.Conv2d(filters[0], 1, 1, 1),
            nn.Sigmoid(),
        )

    def forward(self, x):
        # Encode
        x1 = self.input_layer(x) + self.input_skip(x)#64,240,240
        x2 = self.residual_conv_1(x1)#128,120,120
        x3 = self.residual_conv_2(x2)#256,60,60
        add_x3=self.add_residual_conv(x3)#512,30,30
        # Bridge
        x4 = self.bridge(add_x3)#1024,15,15
        # Decode
        x4 = self.add_upsample(x4)#1024,30,30
        add_x5=torch.cat([x4,add_x3],dim=1)#1024+512,30,30
        add_x6=self.add_up_residual(add_x5)#512,30,30
        add_x6=self.upsample_1(add_x6)#512,60,60
        x5 = torch.cat([add_x6, x3], dim=1)#768,60,60

        x6 = self.up_residual_conv1(x5)#256,60,60

        x6 = self.upsample_2(x6)#256,120,120
        x7 = torch.cat([x6, x2], dim=1)#384,120,120

        x8 = self.up_residual_conv2(x7)#128,120,120

        x8 = self.upsample_3(x8)#128.240.240
        x9 = torch.cat([x8, x1], dim=1)#192,240,240

        x10 = self.up_residual_conv3(x9)#64,240,240

        output = self.output_layer(x10)#1,240,240

        return output
if __name__ == '__main__':
    model= ResUnet(1)
    print(model)