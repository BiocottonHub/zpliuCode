# :baby_chick: 可重复使用脚本合集



+ [:currency_exchange: log2求表达差异](#log2awk)



***



#### log2.awk

使用方法：

+ `fild1`指定A基因表达量所在列
+ `fild2`指定B基因表达量所在列
+ 最终会在最后一列增加，表达倍数的log2值

```bash
awk -v fild1=7 -v fild2=8 -f ~/github/zpliuCode/script/log2.awk  表达量文件 >输出文件
```

