from random import gauss
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
import math


class GreedyOutlier:
    number_of_tuple = 1000
    num_of_attr = 5
    num_of_InEq = 2 * number_of_tuple
    selected_tuple_index = 0
    selected_attr = 3
    membership_threshold = 0.5

    # set_attr = [0, 1, 3]
    # alpha_attr = [1, 1, 1]
    set_attr = [3]
    alpha_attr = [1]


    def __init__(self,table = None):
        if table == None:
            self.create_table()
        print self.table
        self.normalize_t()
        self.value()
        val = []
        for i in self.table:
            val.append(i[self.set_attr[0]])
        print val

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

    def current_outlier_score(self, list_of_value):
        outlier_score=0
        for selected_index in range(len(list_of_value)):
            outlier_score += (self.tuple_val_m[selected_index] - np.mean(list_of_value[selected_index]))**2
        outlier_score = math.sqrt(outlier_score)
        return outlier_score

    def outlier_s(self):
        list_of_value = []
        for attr in self.set_attr:
            tmp = []
            for index in range(self.number_of_tuple):
                if self.inclusion[index] == 1:
                    tmp.append(self.table[index][attr])
            list_of_value.append(tmp)
        outlier_score=0
        for selected_index in range(len(list_of_value)):
            outlier_score += (self.tuple_val_m[selected_index] - np.mean(list_of_value[selected_index]))**2
        outlier_score = math.sqrt(outlier_score)
        return outlier_score

    def if_add(self, list_of_value, value):
        tmp_list_of_value=[]
        for i in range(len(list_of_value)):
            tmp_list_of_value.append(list_of_value[i]+[value[i]])
        return self.current_outlier_score(tmp_list_of_value)

    def greedy_strategy(self):
        pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
        distance_from_given_tuple = []
        self.inclusion = [0] * self.number_of_tuple
        self.inclusion[self.selected_tuple_index] = 1
        for row in range(self.number_of_tuple):
            tmp_dist = 0
            for attr in self.set_attr:
                tmp_dist += (self.tuple_val_m[self.set_attr.index(attr)] - self.table[row][attr]) ** 2
            distance_from_given_tuple.append(math.sqrt(tmp_dist))
        farest_point = distance_from_given_tuple.index(max(distance_from_given_tuple))
        try:
            self.inclusion[farest_point] = 1
        except:
            for far_point in farest_point:
                self.inclusion[far_point] = 1
        number_of_try = 1
        while(self.is_more_tuple() != -1):
            # print "This is the inclusion vector: "
            # print self.inclusion
            # print "This is outlier score for this set: "
            o_score = self.outlier_s()
            # print o_score
            self.histogram(pdf, number_of_try, o_score)
            number_of_try += 1
        pdf.close()

    def is_more_tuple(self):
        isThere = -1
        list_of_value = []
        for attr in self.set_attr:
            tmp = []
            for index in range(self.number_of_tuple):
                if self.inclusion[index] == 1:
                    tmp.append(self.table[index][attr])
            list_of_value.append(tmp)
        current_score = self.current_outlier_score(list_of_value)
        index_of_candidates = []
        diff_value_of_candidates = []
        for index in range(self.number_of_tuple):
            if self.inclusion[index] == 0:
                this_index_val = []
                for attr in self.set_attr:
                    this_index_val.append(self.table[index][attr])
                if_add_this_index_score = self.if_add(list_of_value, this_index_val)
                if if_add_this_index_score > current_score:
                    index_of_candidates.append(index)
                    diff_value_of_candidates.append(if_add_this_index_score - current_score)
        if len(index_of_candidates) > 0:
            max_index = diff_value_of_candidates.index(max(diff_value_of_candidates))
            self.inclusion[index_of_candidates[max_index]] = 1
            isThere = 1
        return isThere

    def histogram(self, pdf, number_of_try, outlier_score):
        for attr in self.set_attr:
            selected_values = []
            not_selected_values = []
            all_value=[]
            for ind in range(self.number_of_tuple):
                all_value.append(self.table[ind][attr])
                if self.inclusion[ind] == 1:
                    selected_values.append(self.table[ind][attr])
                else:
                    not_selected_values.append(self.table[ind][attr])
            fig = plt.figure(ind)
            y_selected_values = [1] * sum(self.inclusion)
            plt.scatter(selected_values, y_selected_values, color='r')

            y_not_selected_values = [0] * len(not_selected_values)
            plt.scatter(not_selected_values, y_not_selected_values, color='b')

            y_all_values = [0.5] * len(all_value)
            plt.scatter(all_value, y_all_values, color='c')

            plt.title("The "+str(number_of_try)+"'th try for attribute "+str(attr)+" scatter plot with outlier score "+str(outlier_score))
            plt.xlabel("Value")
            plt.ylabel("Frequency")
            plt.axvline(x=self.tuple_val_m[self.set_attr.index(attr)], color='k')
            pdf.savefig(fig)
            plt.draw()
            plt.pause(0.01)



a = GreedyOutlier()
a.greedy_strategy()
