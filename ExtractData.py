# 打开文件，第一个参数为文件的路径，第二个参数为读写权限，r代表只允许读
f = open('data/train.tag', 'r')

# 读取文件的每一行
for line in f.readlines():

    # 若该行包含“_”，则证明该行不是序号，是数据
    if "_" in line:
        # print(line)

        # 将每一行按照空格分隔，则会返回一个列表，列表内的元素为 单词_标签 的形式
        line = line.split(" ")
        print(line)

        for word_tag in line:
            word_tag = word_tag.split("_")
            word = word_tag[0]
            tag = word_tag[1]
            print(word)
            print(tag)
            # 后面可以加一些对word和tag的操作
