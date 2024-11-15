class Box:
    """ Handles box operations"""
    def __init__(self, Arg, centering=(0, 0, 0)):

        self.box = self._MakeBox(Arg, centering)


    @classmethod
    def FromLoEndAndSize(cls, LoEnd, Size, centering=(0, 0, 0)):
        Arg = list(LoEnd)
        for Lo, S in zip(LoEnd, Size):
            Arg.append(Lo + S - 1)
        Arg.append('Box')
        return cls(Arg, centering)

    @classmethod
    def FromAnotherBox(cls, box):
        return cls(box.LoEnd() + box.HiEnd()+('Box',), centering=box.Centering())





    def ContainsSlice(self, d, dir):
        return self._BoxContainsSlice(self.box, d, dir)

    def ContainsPoint(self,s):
        return self._BoxContainsPoint(self.box, s)

    def ContainsBox(self,b1):
        return self._BoxContainsBox(self, b1)

    def HiEnd(self):
        return self.box[1]

    def LoEnd(self):
        return self.box[0]

    def Centering(self):
        return self.box[2]

    def Size(self):
        return [max(0, h - l + 1) for l, h in zip(self.LoEnd(), self.HiEnd())]

    def Ncells(self):
        n = 1
        for s in self.Size():
            n *= s
        return n

    def Grow(self, Ghosts):
        box=list()
        box.append(tuple([l - g for g,l in zip(Ghosts,self.LoEnd())]))
        box.append(tuple([l + g for g,l in zip(Ghosts,self.HiEnd())]))
        self.box = (box[0], box[1], self.box[2])

    def IsEmpty(self):
        return any(x == 0 for x in self.Size())

    def Intersection(self, box):
        if self.Centering() != box.Centering():
            raise ValueError("Boxes with different centering cannot be intersected")
        Lo = []
        Hi = []
        if self.IsEmpty() or box.IsEmpty():
            Lo = [0 for x in self.LoEnd()]
            Hi = [-1 for x in self.HiEnd()]

            Args=[]
            for l in Lo:
                Args.append(l)
            for h in Hi:
                Args.append(h)

            return Box(Args)


        for l1, l2, h1, h2 in zip(self.LoEnd(), box.LoEnd(), self.HiEnd(), box.HiEnd()):
            if l1 <= l2 <= h2 <= h1: #segment ob box is contained in self
                Lo.append(l2)
                Hi.append(h2)
                continue
            if l2 <= l1 <= h1 <= h2: # segment of self is contained in box
                Lo.append(l1)
                Hi.append(h1)
                continue
            if h1 < l2 or h2 < l1: # emtpy intersection
                Lo.append(0)
                Hi.append(-1)
                continue
            if l1 <= l2 <= h1 <= h2:
                Lo.append(l2)
                Hi.append(h1)
                continue
            if l2 <= l1 <= h2 <= h1:
                Lo.append(l1)
                Hi.append(h2)
                continue



        Args=[]
        for l in Lo:
            Args.append(l)
        for h in Hi:
            Args.append(h)
        Args.append('Box')
        return Box(Args, centering=self.Centering())






    def _MakeBox(self,Args, centering=(0,0,0)):
        """ make a box as a tuple of low_end, high_end"""
        if Args[-1] != 'Box':
            raise ValueError(" Box constructor called with non box argument")

        low_end = []
        high_end = []

        for i in range(len(Args[:-1]) // 2):
            low_end.append(Args[i])
            high_end.append(Args[i + len(Args[:-1]) // 2])

        return (tuple(low_end), tuple(high_end), tuple(centering[0:len(low_end)]))

    def _BoxContainsSlice(self, b, s, dir):
        """ returns true if b.lo[dir]<=s<=b.hi[dir]"""
        return (b[0][dir] <= s and s <= b[1][dir])

    def _BoxContainsPoint(self, b, s):
        """ returns true if s is in box """
        out=True
        for (d, index) in enumerate(s):
            out &= self.ContainsSlice(index, d)
        return out

    def _BoxContainsBox(self, b0, b1):
        """ returns true if b1 is contained in b"""
        r = True

        for (l0, h0,l1,h1) in zip(b0.box[0], b0.box[1], b1.box[0], b1.box[1]):
            r = r and (l0 <= l1 and h0 >= h1)
        return r

    def __str__(self):
        return "Low End %s, High End %s, Centering %s " % (self.box[0], self.box[1], self.box[2])

    def toCPP(self):
        return (tuple(self.LoEnd()+self.HiEnd()),)