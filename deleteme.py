
class polygon:
    def __init__(self, side):
        self.side = side

class rectangle(polygon):
    def __init__(self, side, otherside=2):
        polygon.__init__(self, side)

class square(rectangle):
    def __init__(self, side, otherside):
        rectangle.__init__(self, side)

sq = square(1, 2)
