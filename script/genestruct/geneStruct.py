from readgtf import getGeneInfo
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.ticker as tick
import sys
'''
刻度上显示所有数字
'''


# def currency(x, pos):
#     return '{:.0f}'.format(x)


geneDict = getGeneInfo(sys.argv[1])
geneId = sys.argv[2]
geneObject = geneDict[geneId]
Num = geneObject.getTranscriptNum()  # 转录本数目
stand = geneObject.getChain()  # 获取正负链
[start, end] = geneObject.getMaxLengthTranscript()  # 获取最长的转录本显示范围

###########################################
'''
开始绘制图片
+ 设置x轴和y轴显示范围
+ 隐藏刻度线及坐标轴

'''
###########################################
interval = 0.5  # 控制转录本间的间隔
fig = plt.figure(1)
ax = fig.add_axes([0, 0.1, 0.85, 0.85])
start = start-200
end = end+200
ax.set_xlim([start, end])
ax.set_ylim([1-interval, Num+interval])

# 设置y轴显示范围以及显示的label
ax.set_yticks(np.arange(1, 1+interval*Num, interval))  # 设置多少个刻度
ax.set_yticklabels([geneObject.getTranscriptName(i) for i in range(Num)])
# 控制坐标轴是否显示, spines 骨脊
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
# 控制刻度以及label显示在哪一边
# ax.xaxis.set_major_formatter(tick.FuncFormatter(currency))  # 实例化formater函数
if stand == "-":
    ax.get_yaxis().tick_left()
    arrow_style = '->'
    ax.set_xticks(np.linspace(start, end, 6))  # 控制刻度数
    print(np.linspace(end, start, 6))
    ax.set_xticklabels(['{:.0f}'.format(i)
                        for i in np.linspace(end, start, 6)])  # 反着输出刻度值
else:
    ax.get_yaxis().tick_left()
    arrow_style = "->"
    # print(start)
    # print(np.linspace(start, end, 6)[1])
    ax.set_xticks(np.linspace(start, end, 6))
    ax.set_xticklabels(['{:.0f}'.format(i)
                        for i in np.linspace(start, end, 6)])


ax.get_xaxis().tick_bottom()
# 隐藏y轴刻度只显示label，x轴则显示刻度并且朝外
ax.tick_params(axis=u'y', which=u'both', length=0)
ax.tick_params(axis=u'x', which=u'both', direction='out')

for i in range(Num):
    transcripCoordinate = geneObject.getTranscriptCoordinate(i)
    if stand == "+":
        arrowLocation = [
            (transcripCoordinate[0]-200, transcripCoordinate[-1]+200), (interval*i+1, interval*i+1)]
        (arrow_start, arrow_end) = zip(*arrowLocation)
    else:
        arrowLocation = [
            (end-transcripCoordinate[-1]+start-200, end-transcripCoordinate[0]+start+200), (interval*i+1, interval*i+1)]
        (arrow_start, arrow_end) = zip(*arrowLocation)
    arrow = mpatches.FancyArrowPatch(
        arrow_start,
        arrow_end,
        arrowstyle=arrow_style,
        mutation_scale=20,  # 箭头缩放比例
        lw=1,
        color='b',
        antialiased=True
    )
    ax.add_patch(arrow)
    # 注释exon元件,也分正负链
    exonCordinate = geneObject.getTranscriptCoordinate(i)
    if stand == "+":
        for index in range(0, len(exonCordinate), 2):
            element = [(exonCordinate[index], interval*i+1),
                       (exonCordinate[index+1], interval*i+1)]
            (line_x, line_y) = zip(*element)
            ax.add_line(
                lines.Line2D(
                    line_x,
                    line_y,
                    linewidth=5,
                    solid_capstyle='butt',
                    solid_joinstyle='miter',
                    antialiased=False,
                    color='red'
                )
            )
    else:
        for index in range(0, len(exonCordinate), 2):

            element = [(start+end-exonCordinate[index], interval*i+1),
                       (start+end-exonCordinate[index+1], interval*i+1)]  # 颠倒一下顺序
            (line_x, line_y) = zip(*element)
            ax.add_line(
                lines.Line2D(
                    line_x,
                    line_y,
                    linewidth=5,
                    solid_capstyle='butt',
                    solid_joinstyle='miter',
                    antialiased=False,
                    color='red'
                )
            )
# 给指定区域填充颜色


def fillRegion(reg_start, reg_end):
    if stand == "+":
        ax.fill([reg_start, reg_start, reg_end, reg_end], [1-interval,
                                                           Num+interval,  Num+interval, 1-interval, ], facecolor='y', alpha=0.2)
    else:
        ax.fill([end-reg_start+start, end-reg_start+start, end-reg_end+start, end-reg_end+start], [1-interval,
                                                                                                   Num+interval,  Num+interval, 1-interval, ], facecolor='y', alpha=0.2)


# fillRegion(29131397,  29131751)

fig.savefig(sys.argv[2]+'.pdf', dpi=150, bbox_inches="tight")
