from cvxopt import matrix,solvers
from random import gauss
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
import math


class MaximumContextOutlier:
    number_of_tuple = 1000
    num_of_attr = 5
    num_of_InEq = 2 * number_of_tuple
    selected_tuple_index = 0
    selected_attr = 3
    membership_threshold = 0.5

    set_attr = [0, 1, 3]
    alpha_attr = [1, 1, 1]
    # set_attr = [3]
    # alpha_attr = [1]


    def __init__(self,table = None):
        if table == None:
            self.create_table()
        # # print self.table
        # self.normalize_t()
        # # print self.table
        self.value()
        self.epsilon_definer()
        self.create_p0()
        self.create_q0()
        self.create_g()
        self.create_h()
        self.create_A()
        self.create_b()


    def normalize_t(self):

        for attribute in self.set_attr:
            column_m = []
            for row in self.table:
                column_m.append(row[attribute])

            min_selected_attr = min(column_m)
            diff_selected_attr = max(column_m) - min_selected_attr
            for r in self.table:
                r[attribute] = (r[attribute] - min_selected_attr) / diff_selected_attr

    def value(self):
        self.tuple_val_m = []
        for a in self.set_attr:
            self.tuple_val_m.append(self.table[self.selected_tuple_index][a])
        self.colmnSum()

    def create_table(self):
        table = []
        for i in range(self.number_of_tuple):
            tmp_row = []
            for attr in range(self.num_of_attr):
                tmp_row.append(gauss(attr, 1))
            table.append(tmp_row)
            # table[self.selected_tuple_index][0] = 4
            # table[self.selected_tuple_index][1] = 6
            table[self.selected_tuple_index][3] = 4
            self.table = table

        return table

    def colmnSum(self):
        pow2col = []
        pow1col = []
        for attribute in self.set_attr:
            sum2 = 0
            sum = 0
            for row in self.table:
                sum += row[attribute]
                sum2 += row[attribute] ** 2

            pow1col.append(sum)
            pow2col.append(sum2)
        self.pow1col = pow1col
        self.pow2col = pow2col
        print pow2col
        print pow1col

    def create_p0(self):

        p = []
        p0_t_attrs = []

        for attr in range(len(self.set_attr)):
            column_m = []
            for row in self.table:
                column_m.append(row[attr])
            p0_t_attrs.append(self.create_p0_t(column_m))


        for row in range(self.number_of_tuple):
            tmp_row = [0] * self.number_of_tuple
            for attr in range(len(self.set_attr)):
                tmp_row = map(sum,zip([x*self.alpha_attr[attr] for x in p0_t_attrs[attr][row]] ,tmp_row))
            # tmp_row = [x * .5 for x in tmp_row]
            p.append(tmp_row)
        self.p = matrix(p)

        return p

    def create_p0_t(self, selected_attribute_column):
        p0_t = []
        middle_vector = []
        for cell in range(self.number_of_tuple):
            middle_vector.append(selected_attribute_column[self.selected_tuple_index]-selected_attribute_column[cell])
        for i in range(self.number_of_tuple):
            tmp = [0] * self.number_of_tuple
            for j in range(self.number_of_tuple):
                tmp[j] = middle_vector[i] * middle_vector[j]
            p0_t.append(tmp)
        return p0_t

    def create_q0(self):
        q = []
        for i in range(self.number_of_tuple):
            q.append(0.0)
        self.q = matrix(q)
        return q

    def create_g(self):
        g = []

        for cnt in range(self.number_of_tuple):
            tmp_row = [0.0] * self.number_of_tuple
            tmp_row[cnt] = 1.0
            g.append(tmp_row)

        for cnt in range(self.number_of_tuple):
            tmp_row = [0.0] * self.number_of_tuple
            tmp_row[cnt] = -1.0
            g.append(tmp_row)
        tmp = matrix(g)
        self.g = matrix.trans(tmp)

        return g

    def create_A(self):

        a = [1.0] * self.number_of_tuple
        tmp = matrix(a)
        self.a=matrix.trans(tmp)
        return a

    def create_b(self):

        b = [math.floor(self.epsilon * self.number_of_tuple)]
        self.b = matrix(b)
        return b

    def create_h(self):

        h1 = [1.0] * self.number_of_tuple
        h1[self.selected_tuple_index] = 1 - self.membership_threshold
        h2 = [0.0] * self.number_of_tuple
        # h3 = [self.epsilon * self.number_of_tuple]
        h = h1 + h2
        self.h = matrix(h)
        return h

    def solv(self):
        return solvers.qp(self.p, self.q , self.g, self.h,self.a,self.b, kktsolver = 'chol')

    def run(self):

        s = self.solv()
        print sum(list(s['x']))
        print self.epsilon
        membership_score = []
        for score in list(s['x']):
            membership_score.append(1 - score)
        self.histogram(membership_score)
        print sum(membership_score)
        return s

    def histogram(self,opt_values):
        pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
        for i in range(len(self.tuple_val_m)):
            attr_m_values = []
            for r in self.table:
                attr_m_values.append(r[self.set_attr[i]])
            fig=plt.figure(i)
            plt.hist(attr_m_values,bins = 20)
            plt.title("Attribute " + str(self.set_attr[i]) + " Histogram")
            plt.xlabel("Value")
            plt.ylabel("Frequency")
            plt.axvline(x=self.tuple_val_m[i], color='k')
            pdf.savefig(fig)
            fig.savefig('hist ' + str(i) + '.png')
        for i in range(len(self.tuple_val_m)):
            fig=plt.figure(i+len(self.tuple_val_m))
            plt.title("Attribute " + str(self.set_attr[i]) + " Membership Function")
            attr_m_values = []
            for r in self.table:
                attr_m_values.append(r[self.set_attr[i]])
            plt.scatter(attr_m_values,opt_values)
            fig.savefig('mem ' + str(i) + '.png')
            pdf.savefig(fig)

        pdf.close()

    def epsilon_definer(self):
        mins = []
        maxs = []
        val = []
        for i in range(len(self.set_attr)):
            attr_m_values = []
            for r in self.table:
                attr_m_values.append(r[self.set_attr[i]])
            val.append(self.table[self.selected_tuple_index][self.set_attr[i]])
            mins.append(min(attr_m_values))
            maxs.append(max(attr_m_values))
        candidate_val = []
        for i in range(len(val)):
            candidate_val.append(min(
                abs(maxs[i] - val[i]) / abs(maxs[i] - mins[i]), abs(mins[i] - val[i]) / abs(maxs[i] - mins[i])
            ))
        mean = np.mean(candidate_val)
        if len(candidate_val)>0:
            var=np.var(candidate_val)
            self.epsilon = abs(gauss(mean, var + 0.05))
            # self.epsilon = max(candidate_val)
        else:
            self.epsilon = max(candidate_val)


a = MaximumContextOutlier()
s = a.run()
