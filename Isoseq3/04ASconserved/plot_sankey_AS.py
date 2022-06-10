'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-01-05 09:09:40
LastEditors: zpliu
LastEditTime: 2021-01-05 16:04:38
@param: 
'''
import plotly.graph_objects as go
import plotly
pyplt = plotly.offline.plot
#import numpy as np

source = [0, 0, 0, 0,  # RI
          1, 1, 1, 1,  # SE
          2,  # conserved RI
          3,  # intron
          4,  # exon
          5,  # other
          6,  # conservedSE
          7,  # exon
          8,  # intron
          9  # other

          ]
target = [10,  # conserved RI
          11,  # intron
          12,  # exon
          13,  # other
          14,  # conserved SE
          15,  # exon
          16,  # intron
          17,  # other
          18, 18, 18, 18,  # RI
          19, 19, 19, 19  # SE

          ]
# link width
value = [3883, 8314, 1252, 3266,  # RI2
         528, 1207, 709, 788,  # SE2
         3725, 5680, 1021, 2079,  # 2RI
         522, 719, 515, 539  # 2SE
         ]

intronColor = '#FEF3C7'  # intron
exonColor = '#A6E3D7'  # exon
IRColor = '#EBBAB5'  # IR
SEColor = '#F7DC6F'  # SE
OEColor = '#CBB4D5'  # other

color_node = [
    IRColor, SEColor,
    IRColor, intronColor, exonColor, OEColor,
    SEColor, exonColor, intronColor, OEColor,
    IRColor, intronColor, exonColor, OEColor,
    SEColor, exonColor, intronColor, OEColor,
    IRColor, SEColor
]
color_link = [
    IRColor, intronColor, exonColor, OEColor,
    SEColor, exonColor, intronColor, OEColor,
    IRColor, intronColor, exonColor, OEColor,
    SEColor, exonColor, intronColor, OEColor
]
label = [
    'IR', 'ES',
    'IR', 'CI', 'CE', 'OE',
    'ES', 'CE', 'CI', 'OE',
    # allpolyploid
    'IR', 'CI', 'CE', 'OE',
    'ES', 'CE', 'CI', 'OE',
    'IR', 'ES'
]
link = dict(source=source, target=target, value=value, color=color_link)
node = dict(
    label=label,
    pad=40,  # 设置节点的垂直距离
    thickness=10,  # 设定节点厚度
    color=color_node,  # 设置节点颜色
    line=dict(color="#303841", width=0.5),  # 设置边框颜色和厚度
    # x=[0.1, 0.1,
    #    0.1, 0.1, 0.1, 0.1, 0.1,
    #    0.1, 0.1, 0.1, 0.1, 0.1,
    #    0.5, 0.5, 0.5, 0.5, 0.5,
    #    0.5, 0.5, 0.5, 0.5, 0.5,
    #    0.1, 0.1, 0.1
    #    ],  # 节点X的位置，不能从0开始


)
data = go.Sankey(link=link, node=node)
fig = go.Figure(data)
fig.update_layout(
    font=dict(size=10, color='#52616b'),
)
fig.write_image("test.svg",
                width=400,  # 宽度
                height=400,
                engine='kaleido'
                )

#pyplt(fig, filename='export-image.html')
