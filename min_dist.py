from cvxopt import matrix,solvers
from random import gauss
import matplotlib.pyplot as plt


class MaximumContextOutlier:

    number_of_tuple = 50
    num_of_attr = 5
    num_of_InEq = 2 * number_of_tuple
    alpha = 0.5
    index_tuple = 0
    selected_attr = 3
    membership_threshold = 0.7
    epsilon = 0.3

    def __init__(self,table = None):
        if table == None:
            self.create_table()
        # print self.table
        self.normalize_t()
        # print self.table

        self.create_p0()
        self.create_q0()
        self.create_g()
        self.create_h()

    def normalize_t(self):

        column_m=[]
        for r in self.table:
            column_m.append(r[self.selected_attr])

        min_selected_attr=min(column_m)
        diff_selected_attr=max(column_m)-min_selected_attr
        for row in self.table:
            row[self.selected_attr]=(row[self.selected_attr]-min_selected_attr)/diff_selected_attr
        self.tuple_val_m = self.table[self.index_tuple][self.selected_attr]

    def create_table(self):
        table = []
        for i in range(self.number_of_tuple):
            tmp_row=[]
            for attr in range(self.num_of_attr):
                tmp_row.append(gauss(attr, 1))
            table.append(tmp_row)
            table[self.index_tuple][self.selected_attr]=4
            self.table = table

        return table

    def create_p0(self):
        p = []
        for row in range(self.number_of_tuple):
            tmp_row=[0]*self.number_of_tuple
            tmp_row[row]=2*(-self.tuple_val_m**2-self.table[row][self.selected_attr]**2+2*self.tuple_val_m*self.table[row][self.selected_attr])
            p.append(tmp_row)
            self.p = matrix(p)
        print p
        return p

    def create_q0(self):
        q=[]
        sum_attr_m=sum(self.table[:][self.selected_attr])
        for i in range(self.number_of_tuple):
            q.append(2*(self.tuple_val_m ** 2)*self.number_of_tuple-2*self.tuple_val_m*self.number_of_tuple*self.table[i][self.selected_attr]-2*self.tuple_val_m*sum_attr_m+2*sum_attr_m*self.table[i][self.selected_attr])
            self.q = matrix(q)
        return q

    def create_g(self):
        g=[]

        for cnt in range(self.number_of_tuple):
            tmp_row=[0.0]*self.number_of_tuple
            tmp_row[cnt]=1.0
            g.append(tmp_row)

        for cnt in range(self.number_of_tuple):
            tmp_row=[0.0]*self.number_of_tuple
            tmp_row[cnt-1]=-1.0
            g.append(tmp_row)
        tmp_row = [1.0] * self.number_of_tuple
        g.append(tmp_row)
        tmp =matrix(g)

        self.g = matrix.trans(tmp)
        return g

    def create_h(self):

        h1 = [1.0]*self.number_of_tuple
        h1[self.index_tuple] = 1-self.membership_threshold
        h2 = [0.0]*self.number_of_tuple
        h3=[(1 - self.epsilon) * self.number_of_tuple]
        h = h1+h2+h3
        self.h = matrix(h)

        return h

    def solv(self):
        return solvers.qp(self.p, self.q, self.g, self.h, None, None)


    def run(self):
        s = self.solv()
        print list(s['x'])
        self.histogram(list(s['x']))
        return s

    def histogram(self,opt_values):
        attr_m_values=[]
        for r in self.table:
            attr_m_values.append(r[self.selected_attr])
        plt.hist(attr_m_values,bins = 20)
        plt.title("Attribute M Histogram")
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        plt.axvline(x=self.tuple_val_m,color='k')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(attr_m_values,opt_values)

        plt.show()

a = MaximumContextOutlier()
s=a.run()