Example:
���룺gcc gen_ped.c �Co gen_ped
������mpicc kmp.c �Co kmp
���У���������gen_ped����ģʽ����gen_ped Strlen Pedlen Seed Pattern_File������Strlen����ģʽ���ĳ��ȣ�Pedlen����ģʽ������С���ڳ��ȣ�Seed���������ʹ�õ���������Pattern_File���������ݴ洢���ļ���������kmp.c�й̶�ָ�����ļ���Ϊpattern.dat��������ʹ�������µĲ�����
	gen_ped 3 2 1 pattern.dat
����֮�����ʹ������ mpirun �Cnp SIZE kmp m n�����иô�ƥ���������SIZE����ʹ�õĴ�����������m��ʾ�ı������ȣ�nΪ�ı��������ڳ��ȡ���ʵ����ʹ����SIZE=3����������m=18��n=3��
����mpirun �Cnp 3 kmp 18 2
���н����
�洢��pattern.dat�е�ģʽ��Ϊ��qmq
�洢��match_result�е�ƥ����Ϊ��
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
˵����������ʵ���У����ı�������Ϊ18������������ı���Ϊasasasqmqmqmypypyp���ֲ���3���ڵ��ϣ�ģʽ������Ϊ3�����������ģʽ��Ϊqmq����󣬽ڵ�1�ϵõ�����ƥ��λ�ã���+��ʾ������

