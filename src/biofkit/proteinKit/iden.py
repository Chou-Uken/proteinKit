#我在干什么 作者是我 现在在写 我是个小白 不知道在干什么


'''
很多注释
'''

class Book():

    title1: str = ''
    author: str = ''
    page: int = 1

    # create a book
    def __init__(self, title: str, author, page: int):
        self.title1 = title #参数里的title和class下的那行的title1没关系，只是个占位符，self.title1意思是class下一行的title。title是参数里的title）
        self.author = author
        self.page = page

    '''#没有的重载：一个函数名用n次但是不同的参数
    def __init__(self, title, author):
        self.tltle1 = title
        self.author = author
        self.page = -1
    #py不需要析构函数'''

    def pagecut(self, times:int) -> str:
        newpage = self.page//times
        return ("GG, a book into " + str(newpage) +" pages.")


if  (__name__ =='__main__'):
    book3 = Book('4567', '米哈游', 114514)
    print(book3.title1)
    print(book3.pagecut(2))
