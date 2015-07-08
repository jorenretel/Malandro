def generalFactory(ClassOfObject, *args):

    return ClassOfObject.__new__(ClassOfObject, *args)
