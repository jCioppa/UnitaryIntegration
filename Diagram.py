from sage.all import *

class Diagram():

    def __init__(self, sh):
        self._shape = sh
        self._length = len(sh)
        self._size = self.size()

    def shape(self):
        return self._shape


    def size(self):
        size = 0
        for part in self._shape:
            size += part
        return size


    def length(self):
        return len(self.shape())


    def fillings(self):
        d = self.size()
        return [T for T in StandardTableaux(d).list() if T.shape() == self.shape()]


    def fillingPairs(self):
        d = self.size()
        return [S for S in CartesianProduct(StandardTableaux(d).list(), StandardTableaux(d).list()).list() if S[0].shape() == S[1].shape() == self.shape()] 



