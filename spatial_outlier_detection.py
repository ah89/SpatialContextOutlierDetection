from cvxopt import matrix,solvers
from random import gauss
import matplotlib.pyplot as plt
import numpy as np


class MaximumContextOutlier:

    number_of_tuple = 100
    num_of_attr = 5
    num_of_InEq = 2 * number_of_tuple
    alpha = 0.5
    index_tuple = 0
    selected_attr = 3
    membership_threshold = 0.5
    epsilon = 0.5
    set_attr=[0,1,3]
    alpha_attr=[1,1,1]

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

        for attribute in self.set_attr:
            column_m=[]
            for row in self.table:
                column_m.append(row[attribute])

            min_selected_attr=min(column_m)
            diff_selected_attr=max(column_m)-min_selected_attr
            for r in self.table:
                r[attribute]=(r[attribute]-min_selected_attr)/diff_selected_attr
        self.tuple_val_m = []
        for a in self.set_attr:
            self.tuple_val_m.append(self.table[self.index_tuple][a])
        self.colmnSum()

    def create_table(self):
        table = []
        for i in range(self.number_of_tuple):
            tmp_row=[]
            for attr in range(self.num_of_attr):
                tmp_row.append(gauss(attr, 1))
            table.append(tmp_row)
            table[self.index_tuple][0]=4
            self.table = table

        return table

    def colmnSum(self):
        pow2col = []
        pow1col = []
        for attribute in self.set_attr:
            sum2=0
            sum=0
            for row in self.table:
                sum += row[attribute]
                sum2 += row[attribute] ** 2

            pow1col.append(sum)
            pow2col.append(sum2)
        self.pow1col=pow1col
        self.pow2col=pow2col
        print pow2col
        print pow1col

    def create_p0(self):

        p = []
        p0_tm_attrs = []
        p0_t_prod_attrs = []
        p0_tm_column_attrs = []

        for attr in range(len(self.set_attr)):
            p0_tm_attrs.append(self.create_p0_tm(self.tuple_val_m[attr]))
            column_m=[]
            for row in self.table:
                column_m.append(row[attr])
            p0_t_prod_attrs.append(self.create_p0_t_prod(column_m))
            p0_tm_column_attrs.append(self.create_p0_tm_column(column_m))

        for row in range(self.number_of_tuple):
            tmp_row = [0] * self.number_of_tuple
            for attr in range(len(self.set_attr)):
                tmp_row = map(sum,zip(map(sum, zip(p0_tm_attrs[attr][row],p0_t_prod_attrs[attr][row],p0_tm_column_attrs[attr][row])),tmp_row))
            tmp_row = [x * .5 for x in tmp_row]
            p.append(tmp_row)
        tmp = matrix(p)
        PP=np.matrix(p)
        print "P Rank"
        print np.all(np.linalg.eigvals(PP+PP.transpose())>=0)

        self.p = matrix.trans(tmp)
        print self.p
        return p

    def create_p0_tm(self,cell_value):
        p0_tm=[]
        for i in range(self.number_of_tuple):
            tmp=[cell_value]*self.number_of_tuple
            p0_tm.append(tmp)
        return  p0_tm

    def create_p0_t_prod(self,vector_of_attr):
        p0_t_prod=[]
        for i in range(self.number_of_tuple):
            tmp=[0]*self.number_of_tuple
            for j in range(self.number_of_tuple):
                tmp[j]=vector_of_attr[i]*vector_of_attr[j]
            p0_t_prod.append(tmp)
        return p0_t_prod

    def create_p0_tm_column(self,vector_of_attr):
        p0_tm_column=[]
        for i in range(self.number_of_tuple):
            tmp=[vector_of_attr[i]]*self.number_of_tuple
            p0_tm_column.append(tmp)
        return  p0_tm_column

    def create_q0(self):
        q=[]
        for i in range(self.number_of_tuple):
            for attr in range(len(self.set_attr)):
                tmp=self.alpha_attr[attr]*(2*(self.tuple_val_m[attr] ** 2)*self.number_of_tuple-2*self.tuple_val_m[attr]*self.number_of_tuple*self.table[i][self.set_attr[attr]]-2*self.tuple_val_m[attr]*self.pow1col[attr] + 2*self.pow1col[attr]*self.table[i][self.set_attr[attr]])
            q.append(tmp)
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
            tmp_row[cnt]=-1.0
            g.append(tmp_row)
        tmp_row = [1.0] * self.number_of_tuple
        g.append(tmp_row)
        tmp =matrix(g)
        self.g = matrix.trans(tmp)
        PP=np.matrix(g)
        print "G rank"
        print np.linalg.matrix_rank(PP)
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
        return solvers.qp(self.p, self.q, self.g, self.h, kktsolver = 'ldl', options = {'kktreg':1e-6})


    def run(self):
        cnt=0
        print self.p[0, 0]
        for i in range(len(self.set_attr)):
            if self.tuple_val_m[i] > self.pow1col[i]+(self.pow1col[i]**2-self.pow2col[i])**0.5 or self.tuple_val_m[i] < self.pow1col[i]-(self.pow1col[i]**2-self.pow2col[i])**0.5:
                cnt+=1
                print self.alpha_attr[i]*2*(self.tuple_val_m[i]**2+self.pow2col[i]-2*self.tuple_val_m[i]*self.pow1col[i])
        if cnt==0:
            print "This problem has solution in all spatial"
            s = self.solv()
            print list(s['x'])
            self.histogram(list(s['x']))
            return s
        elif cnt>0 and self.p[0,0]>0.0:
            print "This problem does not have solution in all spatial"
            s = self.solv()
            print list(s['x'])
            self.histogram(list(s['x']))
            return s
        else:
            print "This problem does not have solution"

        s = self.solv()
        print list(s['x'])
        self.histogram(list(s['x']))
        return s


    def histogram(self,opt_values):

        # plt.hist(attr_m_values,bins = 20)
        # plt.title("Attribute M Histogram")
        # plt.xlabel("Value")
        # plt.ylabel("Frequency")
        # plt.axvline(x=self.tuple_val_m[],color='k')
        lf=[]

        for i in range(len(self.tuple_val_m)):
            attr_m_values = []
            for r in self.table:
                attr_m_values.append(r[self.set_attr[i]])
            plt.figure(i)
            plt.hist(attr_m_values,bins = 20)
            plt.title("Attribute "+ str(self.set_attr[i]) +" Histogram")
            plt.xlabel("Value")
            plt.ylabel("Frequency")
            plt.axvline(x=self.tuple_val_m[i], color='k')
        for i in range(len(self.tuple_val_m)):
            plt.figure(i+len(self.tuple_val_m))
            attr_m_values=[]
            for r in self.table:
                attr_m_values.append(r[self.set_attr[i]])
            plt.scatter(attr_m_values,opt_values)

        plt.show()

a = MaximumContextOutlier()
s=a.run()