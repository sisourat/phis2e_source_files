Version0: 
two-e-int: 最原始的版本，直接一一对应最初的tptp排列，并且对应qp2的文件排列读取。


Version1:
two-e-int: 在0的基础上修正了逆序ijkl-的所有排序，修改了部分交换，输出，未更改文件名。


Version2:
two-e-int: 在1的基础上进一步修改tptp块，交换了求解，使得GTO1=GTO2和GTO3=GTO4。暂时没有修改文件名。

Version3:
two-e-int: 交换了tptp等的文件名字顺序，依照junwen-SCAOCC的collint_2e逆序排列了相关顺序，并且修改了一些共轭项和相等项的
ijkl的顺序，并且对应qp2的文件排列读取。


Version4:
two-e-int: 在2的基础上，保留了修改后的TPTP块，向下寻找可能的错误。pppp块是否为0也在影响r12项的大小.
设置pppp与tttt块完全相同的结构（仿照junwen）。
QP2中暂时不更改读取文件的ijkl顺序。

Version5:
two-e-int: 在4的基础上，仅仅恢复tptp块到原始的。

Version6:
two-e-int: 在0的基础上，恢复所有，对tttt pppp采取相同排列。

