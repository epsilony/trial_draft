'''
Created on 2012-12-17

@author: epsilonyuan+mb@gmail.com

'''

import numpy as np

class Polygon(object):
    def __init__(self, vertex_xys):
        self.vertex_xys = np.array(vertex_xys,dtype=np.double)
        if len(self.vertex_xys.shape) < 2 or self.vertex_xys.shape[1] < 2 or self.vertex_xys.shape[0] < 3:
            raise ValueError("vertex_xys must have a shape which is bigger than (3,2)")
    
    def ray_crossing(self, xy):
        # Originate from:
        #       Joseph O'Rourke, Computational Geometry in C,2ed. Page 244, Code 7.13
        rcross = 0
        lcross = 0
        n = self.vertex_xys.shape[0]
        xys = self.vertex_xys
        x = xy[0]
        y = xy[1]
        for i in xrange(n):
            if np.alltrue(xys[i] == xy):
                return 'v'
            i1 = (i + n - 1) % n
            
            x_i,y_i = xys[i]
            x_i1,y_i1 = xys[i1]
            rstrad = (y_i > y) != (y_i1 > y)
            lstrad = (y_i < y) != (y_i1 < y)
            
            if rstrad or lstrad:
                xcross = (x_i*y-x_i*y_i1-x_i1*y+x_i1*y_i)/(y_i-y_i1)
                if(rstrad and xcross > x):
                    rcross += 1
                if(lstrad and xcross < x):
                    lcross += 1
        if rcross % 2 != lcross % 2 :
            return 'e'
        if rcross % 2 == 1:
            return 'i'
        else:
            return 'o'
        
def sample_polygon():
    # An triangle
    return Polygon([[0, 0], [1, 0], [0, 1]])

if __name__ == '__main__':
    pg = sample_polygon()
    
