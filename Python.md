1.****Basics**

**Command****：**

​       \# 但行注释，””” 或 ''' 是批量注释

help 

​      ?/??(显示函数源码)

​      help(…)

​      type(…)

​      %history

​      dir(math): 查看函数

​      查看所有内建函数：dir(__builtins__)

查看内建函数的帮助文档：print(print.__doc__)

查看内建模块：

​      import sys >>> sys.builtin_module_names

**Typing short cut**

​       control a/b/f/e

**Check dimensions**    

​      rt.shape

​      data.shape(a)

%run 

%road

%paste

Magic methods:

​    _new_, _init_, _delete_

check file, functions

​     !cat ex1.csv

​     % reset - f      deleted variables

Clean

​      %reset -f 

​      %clear

改变工作环境

import os

os.chdir("/Users/liu/Documents/01Data/01TCGA/01KIRC")

os.getcwd()

path="/Users/liu/Documents/01Data/01TCGA/01KIRC"

os.chdir(path)

or 

Read and save file

​           Read:

df = pd.read_csv(‘ex1.csv’, sep = ‘,’,header = None, names= […], index_col = ‘…’,skiprows = […,…,…])

​              df = pd.read_table('examples/ex3.txt', sep='\s+')

​            Write:

​               df.to_csv(‘out.csv’, index = False, header = False)

​               Filtering out missing data:

​               data.dropna()

​               data[data.notnull()]

​               data.dropna(how=‘all') #only drop rows that are all NA 

df = pd.read_csv('pandas_dataframe_importing_csv/example.csv', names=['UID', 'First   Name', 'Last Name', 'Age', 'Pre-Test Score', 'Post-Test Score’])



Print 

​      data 

​      print(data) 可读性降低。

**Modules**

import numpy as np

import matplotlib.pyplot as plt

import pandas as pd

import seaborn as sans

import statsmodels as sm

  打开文件运行某种操作然后关闭：(实例中是读取文件)

​    with open('/Users/duwangdan/Desktop/test.txt') as file:

​    data = file.read()



​    设置缺失值：

​     sentinels = {'Last Name': ['.', 'NA'], 'Pre-Test Score': [‘.']}

df = pd.read_csv('pandas_dataframe_importing_csv/example.csv', na_values=sentinels)

df

更改环境变量

path = "/Users/liu/Documents/01Data/01TCGA/01BLCA/Cluster"

os.chdir(path)

os.getcwd()

**2****.****data process**

**see:**

​    df.head(10)

​    df1.iloc[0:10,0:10]

​    reviews.loc[:5,”survival”,”status"]

column name:

​     list(data)

list = [‘a string’,’b string’,’c string’]

tuple = (1, 2, 3)

list.append(obj)末尾添加新的对象，

list.extend(seq)在末尾追加一个序列的多个值，

list.insert(1,’a’)索引1处添加a

join:

​     along rows: 

​             df_new = pd.concat([df_a, df_b])

​     along columns: 

​              pd.concat([df_a, df_b], axis=1)



结合两个数据框：

​       pd.merge(df1, df2)

​          <https://jakevdp.github.io/PythonDataScienceHandbook/03.07-merge-and-join.html>

​          <https://chrisalbon.com/python/data_wrangling/pandas_join_merge_dataframe/>

​          Merge two dataframes along the columns:

​       pd.merge(df_new, df_n, on=‘subject_id') # on -> ‘inner’, ‘outer’. how ->         ‘right’, ‘left’

data[]

data.iloc[]

data.loc[]

**3****.****String manipulation**

replace:

rnaSeq1 = re.sub('T','U',dnaSeq)

rnaSeq2 = dnaSeq.replace(’T’,’U’)

Remove:

​      strip(‘…’), rstrip(‘…’), lstrip(‘…’)

Pop:

Insert:

glob

​    glob:* 匹配0或多个字符，？匹配单个字符，[]匹配指定范围的字符

​     glob.iglob:逐个获得匹配的文件



list(‘hello, world’)：把字符串转换成列表

tuple(‘hello,world’):把字符串改成元组

enumerate():输出输出列表与相应的索引，如(0,a) (1,b)…

zip():删除相应的索引，留下元组形势的数据

reversed():反转列表

len():

sorted():

max():

min():

sum():

str():





​                            4.**Visualization**

**plot(x,y1, x, y2, )**

**savefig(‘…’)**

​                            **5. M****athematic algorithm**



logarithms:    import math

​                     math.log(100,10)

display the IPython quick reference card :

​      %quickref 

文件存储

​      df.to_csv(‘example.csv’)

空函数：def nop(): pass, 其中pass语句什么都不做，作用是占位，将来想到了再回来写。

**assignment operators:**常见赋值运算符有 = , += , -= , *= , /=浮点除 , %= 余数, //= 整除, **=指数 

assignment:在python中用一个变量给另一个变量赋值，其实就是给当前内存中的对象添加标签。















lambda函数：g = lambda x,y: x*y

**variables:** 分为global variables 和local variables. 用global强调全局变量。

异常：BaseException, exception, attributeError, IndexError, IOError输入输出, keyboradInterrupt用户中断执行, keyError, NameError, SyntaxError, TypeError, ValueError, ZeroDivisionError

捕捉异常：常用try>except>else>finally模块捕捉异常。如果try语句有异常就执行except除非语句, 如果想捕捉多个异常用多个except句子或者用except(… , … ，…)模块。若无异常执行else语句，无论是否发生异常，finally语句都要执行。

避免关闭文件报错：如果try 语句是f=...，finally语句是f.close() ,那么如果try语句错误，finally语句也不会执行，此时要用 with open(‘…’) as f 语句关闭文件。

字符串的表示：’str’普通字符串，”str”当字符串内用了单引号时，’’’ str ‘’’需要换行的字符串。

转义字符：常用转义字符\0,\a,\b,\t,\n,\v,\f,\r,\e,\”,\’,\\,\(在行尾)，分别是空格、响铃、退格、横向制表符、换行、纵向制表符、换页、回车、转义、双引号、单引号、反斜杠、续航符

raw字符串：r’str’  使转义字符失效。

**运算符：**

标准类型运算符：值比较<,>,<=,>=,==,!=，对象身份比较is,is not，布尔运算not,and,or。

序列类型运算符：1.获取s[i],s[:j], s[i:j], s[i:j:k]从i开始到j之前结束(不包括君j)k步长，2.重复s*n, n*s，3.连接s+t，4.判断x in s, x not in s。

字符串方法：

字符串查找：

find():返回某字符载字符串中第一次出现的位置，如果没有匹配则返回-1。find(‘a’,0,3)中0代表起始位置，3是终止位置。

rfind():返回字符**最后一次**出现的位置，如果没有匹配项则返回-1。

index(): 跟find()一样

rindex()：跟rfind()一样

split():分隔符，taobao.split(‘.’,2) #将某字符串从左开始按照.符分割2次成新的列表，rsplit()从右往左。

统计字符的出现次数：count(), count(‘a’,0,10) #0和10表明起始终止位置。

partition():将字符串从左按照某符号分割成元组。rpartition从右

splitlines(): 按行分割

join():’ ’.join(list1)# 将列表中的字符串按.号拼接。’:’join(str)将字符串的字符按照：号拼接。

strip():str1.strip(‘!’)#去掉strip！两头的！lstrip去左边，rstrip去右边

替换：

love.replace(‘ha’,’la’,3)#将love中三个ha换成la

table= str.maketrans(intab, outtab)_Str1.translate(table)  #intab换成outtab

import re_hii=re.compile(‘Hi’) _  hii.sub(‘hey’,hi) #hi 里的Hi 换成hii。

比较：cmp()将两个字符先转化为ASCII值，然后比较大小，左边比右边小返回-1，左边比右边大返回1，左右相等返回1，该函数仅限python2.  用==比较两个对象的值，用is比较两个对象的id。

大小写：upper()把字符串中的小写字母转为大写字母,isupper()判断字符是否都是大写,lower(),islower(),swapcase()大小写呼唤,capitalize()将字符串首字母大写,title()所有单词首字母大写,istitle()

前后缀：startswith()判断指定符号串是否以指定字符开头,endswith()

判断字符串：isalnum()判断是否由数字和字母组成,isalpha(),isdigit()是否由数字组成,isspace()

对齐与填充：center(10，’+’)居中两侧用+填充至十列,ljust()做对齐,rjust()右duiqi,zfill()右对齐，左侧用0填充至16列

line[]: line[:-1]就是去除文本最后一个字符，返回剩下的部分。

**3****.** **复合语句**

**基本语句：有循环，递归****...if, while, for, try, with****等**

**Loop****循环****:**动静如一，有去无回。动静如一指有么没有变化，要么有同样的变化。while 循环for 循环

for 循环：语法for item in iterating object:  do something。迭代对象可以是string, list, tuple, dictionary, file. 迭代指同一个变量，用不同的数值代替。在序列穷尽时停止。

while循环：while condition：do something。for 和while两者的相同点在于都能循环做一件重复的事情，while循环是在条件不成立时停止。

break, continue, else：在while和for循环中常用语句，将

break是立即终止当前循环，转而执行循环外语句。

continue :跳出循环，进入下一次循环，在while循环中判断条件是否满足，在for循环中判断迭代是否已经结束。

else：如果代码从break处终止，跳出循环，不执行else中的代码，若正常结束循环，则执行else 中的代码。

nested loop:即连着写多个嵌套循环,for row in matrix:  for n in row:….

**Recursion****递归****:**递去-归来，中文比英文要达意，即静中有动，有去有回。基本思想是把规模大的问题转化为规模小的相似的子问题来解决(即数学里的归纳法)。即直接或间接调用自己的自调用函数。递归的代码更简洁，更符合自然逻辑。递归必须有边界条件，比如，n==o 或n==1等等。

**if****：**条件表达式里可以有比较运算符，成员运算符，逻辑运算符。与elif expression2: 和else: 连用，elif 和else不缩进，print 要缩进四个空格。

**try-except-else-finally:**语句可以捕捉异常 ，如果被检测的语句块没有异常，则忽略except后面的语句，多个except子句捕捉多个异常。若无异常就执行else子句。无论是否异常都要执行finally语句。

**with:** 用with语句打开文件可以自动关闭，无需使用f.close()

**input():**键盘输入，输入对象为string格式，返回值也是string格式

**int():**将对象转化为整数格式

**Print():**按照定义的格式输出/打印。常见格式有：%x,%d,%o,%f,%s

print (name1, name2, num1, num2, sep=' -- ', end=‘…n’)  #用- -分离，…，换行。

print ('There are %d punctuation marks.' %(count)) # 输出整数

print ('Numbers are %d and %d' %(num1,num2))

print ('That is %10.2s' %(name2))   # 结果是That is         Li   

print (‘Age:{0:<5d},Height:{1:5.2f} m’) m’.format (age. height)) #Age:26   , Height: 1.61 m

print x, y   #函数可以返回多个值。

str():返回字符串，出来的值是给人看的。

repr()：返回标准字符串，标准字符串就是它在python中的以某种形式被存放，可以通过内建函数eval()重新得到字符串。



**3.Array****数组**

**Three array in Python****：****[list], (tuple),{dictionary}**

**List****列表：**可变，可由不同类型组成   

​       [1,2,’a’,2.1,’君’]也可由元组组成，通过索引进行查找。

Tuple 元组：很随意！！！
​     a = ('foo', 'bar') * 4

​       b = (4, None, 'foo') + (6, 0) + (‘bar’,)





append, insert, remove, pop

列表函数: 函数不会真正改变原列表，只是产生个副本。有如下列表函数：list.count(obj)计算某个元素在列表中出现的次数，list.index(obj)索引位置,pop(obj=list[-1])移除某一元素,,reverse()反向列表中元素，sort([func])改变列表地排序，sorted()不改变列表地排序，enumerate(),copy(),

插入元素：list.append(obj)末尾添加新的对象，如果添加多个元素则已元组的形势存在新列表中，list.extend(seq)在末尾追加一个序列的多个值，

list.insert(1,’a’)索引1处添加a

删除元素：list.pop(-1)跳出并删除列表中索引为-1的元素，默认的话是最后一个元素。

del list[-1]  删除索引为-1的元素，remove(obj)移除某个值的第一匹配项

排序元素：list1.sort(),  sort(reverse=True)从大到小降序，list.reverse()反转, 利用key函数对根据任意元素进行ASCII码升序排列，如Fruit.sort(key=len)#按长度排序, Fruit.sort(key=lambda x: x[3]#按每个元素的第四个字母的ASCII码生序排列。

列表的切片操作：List[1][2] # 列表1中第二个元素。

List comprehensions列表解析:通过设定条件将可迭代的对象解析／转换成一个新的列表。[expression for iter_var in iterable1 …] ,  例如

\>>> [x**2 for x in range(10)] 

[0,1,4,9,16,25,36,49,64,81] 

\>>> [(x+1,y+1) for x in range(2) for y in range(2)]

[(1,1),(1,2),(2,1),(2,2)]

生成器表达式：类似列表解析，但不创建列表，只创建生成器，通过for 循环可以查看生成器的内容。表达式：(expression for iter_var in iterable1….)

range(len()): 全长

list.append(): 变成列表

T**uple****：不可变**

**tuple = (3,4,(4,5))**

比方说只能有sorted(…)函数，不能用sort()函数， [(),(),(),()] 每个元括号里的通常称为列表。如果在程序中以列表的形势传递一个对象的集合，它可能在任何地方改变，如果是用元组的话，则提供了一种完整的约束。通过索引进行查找。

zip(): 将列表变成元组

元组主要用在什么地方？

1.在映射类型中当作键使用。

\>>>def foo(args,args = ‘World!’):   #形式参数

\>>>    print(args1,args2)

\>>>foo(‘Hello,’)  #设置位置参数来开始写实参

Hello,World!

\>>>foo(‘Hello,’, args2=’Python!’)  #设位置参数与关键字参数

Hello,Python!

\>>> foo(args2=’Apple!’, args1=’Hello,’)   #设关键字参数

Hello,Apple!

2.作为函数的特殊类型参数

\>>>def foo(args1,*args):   #*有搜集参数的作用

\>>>    print(args1)

\>>>    print(argst)

\>>>foo(‘Hello’,’Wangdachui’,’Nuyyun’,’Linling’)  #实参，Hello 位置参数，后三个组成元组形成参数

Hello,

(’Wangdachui’,’Nuyyun’,’Linling’) #

3.作为函数的特殊返回值

\>>>def foo()

\>>>    return1,2,3

\>>>foo()

(1,2,3)    #其它语言语言可能只返回1，可是pyhton用元组形式返回所有值，保持完整性。

**Dictionary****字典：**由key和value组合形成映射关系. 键可以是数字，字符串，元组组成。值是随机存储的没有顺序，通过键值进行查找。

创建字典:1.直接创建{‘a ’:0.1,’b’:0.2,’c’:0.3}. 2.创建空字典Dict()={}

dict(): 将列表元组等转转成字典，只要映射关系清楚，系统能都能处理。

{}.fromkeys(seq, val): 用于以seq中元素做字典的键，val为所有值，默认为None。

例题，创建员工信息表时如何将所有员工的工资默认值设置为3000？

\>>>adict={}.fromkeys((‘Wangdachui’,’Niuyun’,’Linling’,’Tianqi’),3000)

dict(zip(…,…)):将元组变成字典。

例题，已知有姓名列表和工资列表，如何生成字典类型的员工信息表？对于有公司代码，公司名称，股票价格的列表，如何构造只含有公司代码与股票价格的字典？

\>>>for i in range(5):   #用循环遍历列表所有元素

\>>>    aStr=pList[i][0]

\>>>    bStr=pList[i][0]

\>>>    aList.append(aStr) 

\>>>    bList.append(bStr)

\>>>aDict=dict(zip(aList,bList))

\>>>print(aDict)

dict[i[0]]=i[2] : 通过键来确定对象，

添加：1.dict[‘…’]=…: 给键赋值。2.Dict1.update(Dict2)把字典2的信息更新到字典1中，无返回值。

dict[‘…’]:查找

‘…’ in FruitDict : 判断是否为字典成员。python2 还可以用has_key().

def func(args1,*argst,**argsd): 字典可以做函数的可变长关键字参数。一个*号是将位置参数收集到一个元组中，两个*号收集到一个字典中。因为函数的长度未定所以是可变长。

\>>>func(‘Hello,’,’Wangdachui’,’Niuyun’,’Linling’,a1=1,a2=2,a3=3)  

字典的内建函数：dict()创建字典,len()键值对个数,hash()哈希函数，返回对象的哈希值。

提取键和值：keys()返回字典中所有key 为列表,values()返回值,get()返回指定键的值,setdefault()跟get一样

删除：1.del dict[‘a’]删除a键与相应值，del dict删除整个字典，2.dict.pop(‘a’)删除a与相应的值并返回a的值. 3.clear()清空并返回空字典

items():遍历字典

copy()：浅复制bDict=aDict.copy()

例子，已知有员工姓名和工资信息表,如何单独输出员工姓名和工资金额？alnfo.keys()和alnfo.values()

例子，人事部门有两份人员和工资信息表，第一份原有信息alnfo，第二份是工资更改和新进人员信息binfo，如何获得完整的信息表? alnfo.update(binfo)

JSON查询：

\>>>x={“name”:”Niuyun”,”address”:

​     {“city”:”Beijing”,”street”:”Chaoyang Road”}

​      }

\>>>x[‘address’][‘street’]

‘Chaoyang Road’

搜索引擎关键词查询：

\>>>import requests

\>>>kw={‘q’:’Python dict’}

\>>>r=requests.get(‘http://cn.bing.com/search’,params=kw)

\>>>r.url

\>>>print(r.text)

**Set****集合**：{‘a’,’b’,’c’,’d’}一个无序不重复的元素的组合. 因此不支持索引，切片等操作。包括set(可变)和frozenset(不可变)两种。

创建集合：set(iterable), frozenset(iterable), set(), frozenset()

| **数学符号**        | **集合运算符** |
| ------------------- | -------------- |
| ∈                   | **in**         |
| ∉                   | not in         |
| **=**               | **==**         |
| **≠**               | **!=**         |
| ⸦**包含于／被包含** | **<**          |
| ⊆**真包含于**       | **<=**         |
| ⊃**包含**           | **>**          |
| ⊇**真包含**         | **>=**         |

| 数学符号                | 关系运算符 |
| ----------------------- | ---------- |
| **∩****交**             | &          |
| ∪并                     | **│**      |
| - 差集                  | -          |
| **∆****对称差集／异或** | ^          |



集合内建函数:

超集子集：a.issubset(b)：返回布尔值,判断a集合是不是b集合子集。a.issuperset(b)：返回布尔值，判断a是不是b的超集。

并交差对称差：union()并,intersection(),difference(),symmetric_difference()

copy()：浅复制

可变集合内建函数: update()把原集合更新成并集,intersection_update()差集,

​                             difference_update(),symmetric_difference_update(),

增加删除集合的元素：add(),remove(),discard(),pop(),clear()



**4.Data acquisition and statistics**

缓冲区：打开文件时，先将文件内容载入缓冲区（缓存），并返回一个指向FILE结构体的指针，接下来对文件的操作，都映射成对缓冲区的操作，只有当强制刷新缓冲区、关闭文件或程序运行结束时，才将缓冲区中的内容更新到文件。就像编辑word文档，并不是立刻将编辑好的内容写入到磁盘上的文件，而是对缓存中的副本进行操作，只有当保存文件时，才将副本同步到磁盘上的文件。

**格式说明符：**b o x c d f e

+m.nf： 整个输出占m列，保留n位数，f表示浮点型

<：左对齐，字符串默认方式，用空格填充右边

\>: 右对齐，整数默认方式

0>5d：用0填充，右对齐，宽度5，

 ^：剧中对其

{{}}：输出{}

转义字符：\0(空字符)，\a(响铃),\b(退格),\t(横向制表符),\n(换行),\v(纵向制表符),\f(换页),\e(转页),\”(双引号),\’(单引号),\\(反斜杠),\(续行符)

dir(str)：用来查看所有字符串方法。

**操作模式**：w,r,wr,rt.   #类Unix平台的换行符是\n，而windows平台用的是\r\n来表示换行，python内部采用的是\n来表示换行符。rt模式下，python在读取文本时会自动把\r\n转换成\n. wt模式下，Python写文件时会用\r\n来表示换行。

**数据获取****:**打开文件才能读写，关闭文件防止系统崩溃。

打开文件：file_obj=open(filename, mode='r', buffering=…)#不写缓冲代表使用默认值-1，0代表不缓冲。

Open 函数的mode参数：r , w, a(文件后加内容),r+,w+,a+.  二进制文件的后边加b，例如rb.

**文件读写：**

file_obj.write (str)

file_obj.read(size) size设计就读几个字符

file_obj.read() 是读之前读剩下的。

除此之外还有：file_obj.readline()，file_obj.readlines()，file_obj.writelines()

file_obj.seek(offset, whence=0):在读写文件的时候都有一个文件指针，数据从文件指针所在的位置开始读写。offset是位移量whence是起始位置。

**文件关闭：****1.f.close(). 2.****一开始就用****with open ('test.txt', 'w') as f**

**文件夹操作：常用到****os****，和****shutill****模块**

常用statement:

查看模块功能：print(dir())

查看当前工作地址：print(os.getcwd())

改变工作环境：os.chdir(‘/Users/coreyschafer/Desktop/‘)

显示目录文件：print(os.listdir())

创建文件夹：os.makedirs(‘folder1/folder2’)  或os.makedir(folder)

删除文件夹/文件：os.rmdir() #只能删除空目录或 os.removedirs(‘’)或os.remove(“file”) 

​                                shutil.rmtree(“dir”)

返回文件名和路径名：basename()   或 dirname()  

返回有路径的文件名：os.path.join(‘d:\\library','book.txt') #'d:\\library\\book.txt'

返回目录名和文件名元组：os.path.split()。#('d:\\library', 'book.txt')

返回文件名和扩展名元组：os.path.splitext()。#('book', ‘.txt')

获取文件名：os.path.basename()

读取和设置环境变量：os.linesep 

重命名：os.rename(old, new)

获取文件属性：os.stat(file)

修改文件权限与时间戳：os.chmod()

获取文件大小：os.path.getsize(filename)

复制文件：shutil.copyfile(“oldfile”，“newfile”)  或者shutill.copy(“oldfile”,”newfile”)

复制文件夹：shutil.copytree(“olddir”,”newdir”)

移动文件：shutil.move(“oldpos”,”newpos”)

文件编辑与中文处理：

**网络数据获取****:**

**数据获取：**本地数据用open,write.closed，网络数据用urllib,requests,scrapy抓取，用beautifulSoup，redex解析

读取CSV文件：有一个NBA数据集，包含运动员和他们在2013-2014赛季的表现，看看如何玩转这些数据。

import pandas as pd   #用pandas库以使用dataframe

nba = pd.read_csv(r'c:\test\axp.csv')  

print(nba.shape)   #查看矩阵的行列数

print(nba,head(1)) #查看第一行



例子，求道指成分股中30只股票最近一次成交价的平均值？股票最近一次成交价大于等于180的公司名？

\>>>djidf.lasttrade.mean()

\>>>djidf[djidf.lasttrade>=180].name

统计美国运通公司近一年股票的涨和跌分别的天数？

\>>>len(quotesdf[quotesdf.close>quotesdf.open])

\>>>len(quotesdf)-123

统计美国运通公司近一年相邻两天收盘价的涨跌情况？

\>>>status=np.sign(np.diff(quptesdf.close))

\>>>status

array([1.,1.,-1,…,-1.,1.,1.])   #1表示涨，-1表示跌

\>>>status[np.where(status==1.)].size

132

\>>>status[np.where(status== -1.)].size

118

例子，按最近一次成交价对30只道指成分股股票进行排序。根据排序结果列出前三甲公司名。

\>>>tempdf=djidf.sort_values(by=’lasttrade’,ascending=false)  #ascending=false就是让其逆序

\>>>tempdf[:3].name #切片操作

统计本年度1月份的股票开盘天数？

\>>>t=quotesdf[(quotesdf.index>=’2017-01-01’) & (quotesdf.index<’2017-02-01’)]  #2017-01-01应该是dataframe的index

\>>>len(t)  #得天数

统计近一年每个月的股票开盘天数？

Import time    #help(time)查看功能

…

Listtemp=[]

For i in range(len(quotesdf)):

​    temp=time.strptime(quotesdf.index[i],”%Y-%m-%d”)

listtemp.append(temp.tm_mon)

tempdf=quotesdf.copy()

tempdf[‘month’]=listtemp

print(tempdf[‘month’].value_counts())

Grouping：分类，汇总

统计近一年每个月的股票开盘天数？

\>>>x=tempdf.groupby(‘month’).count()   #(‘month’)是分组依据

\>>>x.close

统计近一年每个月的总成交量？两种方法

\>>>tempdf.groupby(‘month’).sum().volume  #先求和，在计算成交量

\>>>tempdf.groupby(‘month’). Volume.sum()   #先计算成交量，再求和，更高效

Merge：append追加concat连接 join (circle类型的连接)

把美国运通公司本年度1月1日至1月5日间的股票交易信息合并到近一年中前两天的股票信息中？

\>>>p=quotesdf[:2] #切前两个

\>>>q=quotesdf[‘2017-01-01’:’2017-01-05’]

\>>>p.append(q)   #p是已有的，q是增加的

将美国运通公司近一年股票数据中的前5个和后5个合并。

\>>>pieces=[tempdf[:5],temdf[len(tempdf)-5:]]    #切前五个和后五个

\>>>pd.concat(pieces)  

如何将两个不同逻辑结构的对象连接？

\>>>piece1=quotesdf[:3]  #这里quotesdf是不含月份的股票数据

\>>>piece2=tempdf[:3]  #tempdf是含月份的股票数据

\>>>pd.concat([piece1,piece2],ignore_index=true)

Concat函数还有其它常用参数：objs, join, keys, names, ignore_index, axis, join_axes, levels, verify_inregrity

将美国运通公司和可口可乐公司近一年中每个月的交易总量表(包含公司代码)与30只道琼斯成分股股票信息合并。

\>>>pd.merge(djidf.drop([‘lasttrade’],axis=1),Akdf,on=’code’.) #对于djidf文件仍掉lasttrade列, axis是啥意思不知道。对于Akdf,基于code列。

merge函数的其它参数：left,right,how,on , left_on, right_on, left_index, right_index, sort, suffixes,copy

Control Structures

Eceptions & Files

​    如果当try后的语句执行时发生异常，python就跳到第一个匹配该异常的except子句并执行，异常处理完毕，控制流就通过整个try语句(除非在处理异常时又引发新的异常)

More Types 

Functional Programming

lambdas能嵌入到其他表达式当中的匿名函数，

map & filter 

generators

decorators

Recursion 

Def factorial(x):

  If x == 1:

Return 1

  Else:

Return x * factorial(x-1)

Print(factorial(5))

Result:>>> 120



Sets

Itertools

Module  

​                                                            **Regular Expressions**



**正则表达式：**正则表达式本身是一种小型的、高度专业化的编程语言，做匹配替换的首选，而在python中，通过内嵌集成re模块，我们可以直接调用来实现正则匹配。正则表达式模式被编译成一系列的字节码，然后由用C编写的匹配引擎执行。

**()****捕捉组，****[]****字符集，**

**正则语法：**

**1.****复杂的正则表达式可以用简单的连接而成，如果字符串****p****与****A****匹配，****q****与****B****匹配的话，那么字符串****pq****也会与****AB****匹配。****2.****正则表达式可以包含特殊字符和普通字符。****’A’,’a’****和****’0’****都是普通字符。****’****标点符号****’****是特殊字符。**





**正则匹配：**

**re.match:****从字符串的开头开始匹配，匹配上一个就停止。**

​                **例如****:re.macth(r’c++’,text)**

**re.search:****从字符串的任意位置开始匹配，匹配上一个就停止。**

**re.findall:****返回所有能匹配上的字符串列表。**

**re.finditer:**

**正则替换：**

**例题，搜索全部的****sbs****，然后替换成****sb’s****，初级正则这样****s/\bsbs/s\b's/g**

高级的会这样s/(?<=\bsb)(?=s\b)/'/g

**正则引擎：**

**RegExr****可视化网站****:  https://regexper.com**

Metacharacters : ^ . $ !



**6. advanced python**

**python****包管理工具**：

distutils: python的标准库，经常使用的setup.py就是基于distutils实现的，为开发者提供方便的打包或安装方式。python setup.py sdist 生成安装包，python setup.py install安装，还可以编写setup.py内容。

setuptools / distribute：是distutils的增强，同过ez_setup.py来安装，使用.egg格式文件。安装后使用easy_install这个工具安装.tgz 和.egg格式的文件。通过easy_install —help 学习。

pip：最流行的包管理工具，通过VCS或浏览器访问地址来安装python包，还支持卸载。例如pip uninstall somepackage.  通过两种方法安装pip：1.python get-pip.py .  2. 通过setup.py 进行安装。 通过pip —help学习。

pip常用命令：

pip install somepackage, pip uninstall somepackage, pip list, pip list —outdated, pip install —upgrade somepackage, pip show —files somepackage, pip install somepackage==1.0.4, pip install ‘somepackage>=1.0.4’, pip freeze> requirements.txt, pip install -r requirements.txt

**Scipy:**第三方库，一个软件生态系统，主要有六个科学核心库.  到<http://scipy.org> 下载安装或使用安装好的Anaconda.

Scipy的数据结构：ndarryN维数组, Series变长字典, DataFrame数据框

**1.NumPy:** 它的思维模式是面向数组。优点如下：1.强大的ndarray对象和ufunc函数2.精巧的函数3.适合线性代数和随机数处理等科学计算4.有效的通用多维数据，可定义任意数据类型5.无缝对接数据库。

ndarray:所有元素的类型必须相同。axis轴：每个线性的数组称为一个轴，axis=0表示沿着第0轴进行操作，即对每一列进行操作；axis=1,沿着第一轴进行操作，表示对每一行进行操作。rank秩：一维数组的秩为1，二维数组的秩为2.

创建数组： 创建函数如下：arrange,

array[]: 创建多维数组，array([(0,1,2,3,),(4,5,6,7),(8,9,10,11)]) 

fromstring(): 将字符串转为一维数组，np.fromstring(‘1 2 3 4’,dtype=int, sep=‘,’)

fromiter():从可迭代对象中读取数据,并转成一维数组，iterable=(x*x for x in range(5)) >>>np.fromiter(iterable,float)

fromfunction():在每个坐标上执行函数表达式，

ones：将数组中所有元素填充为1

ones-like:将多维参数转成数组

zeros, zeros-like：同上

empty，empty-like：返回没有初始化内存的数组

full():用指定的值填充np.full((2,3),4) #用4填充

full_like(): np.full_like(a, 2)# 将a用2填充

eye():返回对角线矩阵，np.eye(3, dtype=float, k=1)  #用k来确定位置

identity():创建单位矩阵 np.identity(2, dtype=int)

copy():浅复制

reshape():通过reshape生成的新数组和原始数组共用一个内存，也就是说，若更改其中数组的元素，另一个数组也将发生改变。两种用法：1.data.reshape(2,4) ,2.np.reshape(data,(2,4))

arange():创建一个一维的等差数列数组np.arange(20).reshape(4,5)

linspace():创建一个高级等差数组。np.linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None)

logspace():对数等差数组

geomspace()：指数等差数组

meshgrid():根据坐标轴向量创建坐标轴矩阵。

mgrid():

ogrid():

fromfile():

r:

**ndarray****的切片与索引：**切片功能与python的list一样，由start，stop，step

一维数组：a[1:6],a[1:6:2],a[:5],a[5:],a[::5]#不设起始终止，只设步长

二维：b[0]行,b[0,:]行,b[:,0]列,b[1,1:] 第2行1开始的行,b[::2,::2]b[3,::3]. 怎么去指定的列？

三维数组：

用数组来索引：在二维数组中索引，A[[0,2],[1,4] 去0行1列的交叉和2行4列的交叉。在三位数组中索引，B[0,1,0],[1,0,1],[4,3,2]，不懂！！

布尔型索引：B>10, B[B>10]筛选出A中大于10的元素 ,A[A% ==0]筛选出A中的偶数

列出行、列、元素：for row in B: >>>print(row)列出B的行，B.T得到转置矩阵，B.flat 生成器，B.flatten() 将B转成一维数组

最大最小址索引：np.argmax(b)整个矩阵的最大值索引，np.argmin(a,1)每一行最小值索引，np.argmax(a,1)每一行最大值索引，np.argmax(a,0)每一列最大值,np.argmin(a,0)每一列最小值

**ndarray****的基本运算：****a//2****除****2****的整数****,a%3****除余数****, np.dot(m.n)****点乘**

**numpy****中的统计学：**

np.sum(a), np.min(a), np.max(a), np.max(a,axis=0)按列求最大值, np.sum(a,axis=1)按行求最和

np.mean(a)所有元素的平均值, a.mean(a,0)每一列的平均值，np.average(a,1,weights=[1,0,3])加权平均值，np.mean(a)中位数,np.nanmean(b)忽略忽略nan(not a number)计算平均值，np.var(b)方差，np.std(b)标准差，np.nanstd(b)

**2.Pandas:**  **处理矩阵数据的强大包。**1.基于SciPy library和NumPy。2.拥有高效的Series和DataFrame数据框。3.强大的可扩展数据操作与分析功能。3.高效处理大数据集的切片等功能。5.提供优化库功能、读写多种文件格式，如csv, HDF5 

用列表创建数据框：df=pd.DataFrame([3,4,7,9])

用ndarray创建数据框：df=pd.DataFrame([np.arange(12),reshape(3,4))

用字典创建数据框：

data={‘a’:[1,2],’b’:[2,3],’c’:[5,6]}>>>df=pd.DataFrame(data,index=[3,4,5)#一维的话index=[0]

newdata={‘a’:{‘L1’:11,’L2’:22},’b’:{‘L1’:33}}#使用了嵌套字典第一层字典的键为columns,第二层字典的键为index，第二层字典的值为数据框的值。

打印：df3.columns, df[‘…’],df.index

赋值：df.[‘a’][‘L3’]=44

填充：df[‘b’]=s, 把S填充到df数据框的b位置，系统会用Nan自动对齐

Series: 就是竖起来的list

将列表转换为series：ser1=pd.series([2,3,4.5])

讲元组转换为series:  ser2=pd.series((2,3,-5))

将字典转换为series: dic={‘a’:1,’b’:2,’c’:3}>>>ser3=pd.series(dic)

查看：ser.index, ser.values, ser[[‘a’,’b’]]

赋值：ser[[‘a’,’b’]]=8, ser*3, np.exp(ser)

取索引中包含特定值的行：ens2syn[ens2syn.index==“ENSG00000227232.5"]

取出包含某列特定值的行：ens2syn[ens2syn['gene_symbol'].isin(['DDX11L1','MIR6859-1'])]

使用正则表达式选取符合要求的行：

​                 ens2syn[ens2syn.index.str.contains(r'ENSG0000022')].head()