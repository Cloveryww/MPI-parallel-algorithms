Example:
编译：gcc gen_ped.c –o gen_ped
　　　mpicc kmp.c –o kmp
运行：首先运行gen_ped生成模式串，gen_ped Strlen Pedlen Seed Pattern_File。其中Strlen代表模式串的长度，Pedlen代表模式串的最小周期长度，Seed是随机函数使用的种子数，Pattern_File是生成数据存储的文件，这里在kmp.c中固定指定的文件名为pattern.dat。本例中使用了如下的参数。
	gen_ped 3 2 1 pattern.dat
　　之后可以使用命令 mpirun –np SIZE kmp m n来运行该串匹配程序，其中SIZE是所使用的处理器个数，m表示文本串长度，n为文本串的周期长度。本实例中使用了SIZE=3个处理器，m=18，n=3。
　　mpirun –np 3 kmp 18 2
运行结果：
存储于pattern.dat中的模式串为：qmq
存储于match_result中的匹配结果为：
The Text on node 0 is asasas .
The Text on node 1 is qmqmqm .
The Text on node 2 is ypypyp .
This is the match result on node 0
(0)  -
(1)  -
(2)  -
(3)  -
This is the match result on node 1
(4)  -
(5)  -
(6)  +
(7)  -
(8)  +
(9)  -
This is the match result on node 2
(10)  -
(11)  -
(12)  -
(13)  -
(14)  -
(15)  -
说明：该运行实例中，令文本串长度为18，随机产生的文本串为asasasqmqmqmypypyp，分布在3个节点上；模式串长度为3，随机产生的模式串为qmq。最后，节点1上得到两个匹配位置，由+表示出来。

