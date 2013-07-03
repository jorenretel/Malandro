def generalFactory(ClassOfObject,*args) :
  
  print 'been here'
  
  return ClassOfObject.__new__(ClassOfObject,*args)