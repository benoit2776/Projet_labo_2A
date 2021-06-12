import csv
from Well import Well
from Marker import Marker
import matplotlib.pyplot as plt
import numpy as np


class Biseau:
    def __init__(self, wells_file, correlation_results=None, pendage=None):
        """
        :param wells_file:  list that contains the wells and the depths of all markers
        :param correlation_results:  correlation files (created with WeCo)
        :param pendage: file that contains pendages in the different facies
        """

        ##### information extraction from those files ######

        # list of wells names
        self.wells = []

        # dictionnary that will contain the number of markers for each well
        self.nb_markers = {}

        # dictionnary that will contain depths of each markers for the different wells
        self.depths = {}

        ### reading of the wells list file
        file = open(wells_file, 'r')
        lines = file.readlines()
        currentWell = 0
        currentSection = 0
        for i in range(len(lines)):
            if i == 0:
                a, b, well_number = lines[i].split()
                self.well_number = int(well_number)

            if "Well" in lines[i] and len(lines[i].split()) == 1:
                self.wells.append(lines[i][0: len(lines[i]) - 1])
                currentWell = lines[i][0: len(lines[i]) - 1]
                self.depths[currentWell] = []

            if (len(lines[i].split()) == 2) and ("depth" in lines[i]):
                a, nb_markers = lines[i].split()
                self.nb_markers[currentWell] = nb_markers
                currentSection = "depth"

            if (len(lines[i].split()) == 2) and ("facies" in lines[i]):
                currentSection = "facies"

            if currentSection == "depth" and len(lines[i].split()) == 1:
                self.depths[currentWell].append(lines[i][0:len(lines[i]) - 1])

        self.wells_List = []
        self.id = {}
        for well in self.wells:
            w = Well(well, self.nb_markers[well], self.depths[well])
            w.markers.append(Marker(well, 0))
            self.wells_List.append(w)
            self.id[well] = 0

        ### Reading of the results file

        with open(correlation_results) as results:
            file_reader = csv.reader(results, delimiter=',')
            for row in file_reader:
                well_name, marker, cost = row[0], row[1], row[2]
                if well_name == self.wells[0]:
                    well = self.wells_List[0]
                    other_well = self.wells_List[1]
                    other_well_name = other_well.well_name
                    if cost not in well.cost:
                        well.cost.append(cost)
                        well.dico_cost[cost] = [self.id[well_name]]
                        m = Marker(well_name, self.id[well_name])
                        if not well.marker_presence(m):
                            well.markers.append(m)
                        self.id[well_name] += 1
                    else:
                        m = Marker(other_well_name, self.id[other_well_name])
                        if not other_well.marker_presence(m):
                            other_well.markers.append(m)
                        well.get_marker(self.id[well_name] - 1).addCorrel(
                            other_well.get_marker(self.id[other_well_name]))
                        other_well.get_marker(self.id[other_well_name]).addCorrel(
                            well.get_marker(self.id[well_name] - 1))
                else:
                    well = self.wells_List[1]
                    other_well = self.wells_List[0]
                    other_well_name = other_well.well_name
                    if cost not in well.cost:
                        well.cost.append(cost)
                        well.dico_cost[cost] = [self.id[well_name]]
                        m = Marker(well_name, self.id[well_name])
                        if not well.marker_presence(m):
                            well.markers.append(m)
                        well.get_marker(self.id[well_name]).addCorrel(
                            other_well.get_marker(self.id[other_well_name] - 1))
                        other_well.get_marker(self.id[other_well_name] - 1).addCorrel(
                            well.get_marker(self.id[well_name]))
                        self.id[well_name] += 1
                    else:
                        well.get_marker(self.id[well_name] - 1).addCorrel(
                            other_well.get_marker(self.id[other_well_name] - 1))
                        other_well.get_marker(self.id[other_well_name] - 1).addCorrel(
                            well.get_marker(self.id[well_name] - 1))

    def biseaux(self, n, alpha1=None, alpha2=None):
        """
        :param alpha1: angle between d3 and the horizontal (cf. sketch)
        :param alpha2: angle between d2 and the horizontal (cf. sketch)
        :return: Position of all the discontinuities
        """
        if alpha1 == None:
            alpha1 = 22
        if alpha2 == None:
            alpha2 = 22
        w1 = self.wells_List[0]
        w2 = self.wells_List[1]

        for m in w1.markers:
            m_coord = {"abscisse": 0, "depth": self.depths['WellB5'][m.id]}
            if len(m.correlations) > 1:
                correl = []
                for correl_m in m.correlations:
                    correl.append(correl_m.id)
                mini, maxi = min(correl), max(correl)
                mini = {"abscisse": 750, "depth": self.depths['WellB11'][mini]}
                maxi = {"abscisse": 750, "depth": self.depths['WellB11'][maxi]}

                ### lines equations
                d = {}

                # for d1
                print(m_coord["depth"])
                print(mini["depth"])
                print(maxi["depth"])
                a = (float(m_coord["depth"]) - float(mini["depth"])) / (-750)
                b = float(m_coord["depth"])
                d["d1"] = [a, b]

                # for d2
                a = - alpha2 / 45
                print(mini["depth"])
                b = float(mini["depth"]) - a * 750
                d["d2"] = [a, b]

                # for d3
                a = alpha1 / 45
                print(maxi["depth"])
                b = float(maxi["depth"]) - a * 750
                d["d3"] = [a, b]

                # for d4
                a = (float(m_coord["depth"]) - float(maxi["depth"])) / (-750)
                b = float(m_coord["depth"])
                d["d4"] = [a, b]

                if d["d3"][0] < d["d4"][0]:
                    print("Please choose a more important alpha angle")
                    return 0

                # Lines intersections

                Pts = {}

                # intersection d2-d3 : Point A
                xa = (d["d3"][1] - d["d2"][1]) / (d["d2"][0] - d["d3"][0])
                ya = d["d2"][0] * xa + d["d2"][1]
                Pts["A"] = [xa, ya]

                # intersection d1-d3 : Point B
                xb = (d["d3"][1] - d["d1"][1]) / (d["d1"][0] - d["d3"][0])
                yb = d["d1"][0] * xb + d["d1"][1]
                Pts["B"] = [xb, yb]

                # Point C : marker m
                xc = 0
                yc = float(m_coord["depth"])
                Pts["C"] = [xc, yc]

                # intersection d2-d4 : Point D
                xd = (d["d4"][1] - d["d2"][1]) / (d["d2"][0] - d["d4"][0])
                yd = d["d2"][0] * xd + d["d2"][1]
                Pts["D"] = [xd, yd]

                print(Pts)

                ### discontinuity position

                ## Abscisse of the discontinuity
                # parameters of the normal law
                x1, x2 = max(xa, xb, xc, xd), min(xa, xb, xc, xd)
                x_mu = x2 + (x1 - x2) / 2
                x_sigma = (x1 - x2) / 6

                ## Depth of the discontinuity
                y1, y2 = max(ya, yb, yc, yd), min(ya, yb, yc, yd)
                y_mu = y2 + (y1 - y2) / 2
                y_sigma = (y1 - y2) / 6

                limits = [x1, x2, y1, y2]
                x_discontinuity, y_discontinuity = self.normal_law(x_mu, x_sigma, y_mu, y_sigma, n, limits, d)

                X, Y = np.array([xa, xb, xd, xc]), np.array([ya, yb, yd, yc])
                plt.plot(X, Y, 'x', color='b')
                plt.plot(x_discontinuity, y_discontinuity, 'o', color='r')
                plt.grid()
                plt.show()

        for m in w2.markers:
            m_coord = {"abscisse":750, "depth":self.depths['WellB11'][m.id]}
            if len(m.correlations) > 1:
                correl = []
                for correl_m in m.correlations:
                    correl.append(correl_m.id)
                mini, maxi = min(correl), max(correl)
                mini = {"abscisse":0, "depth":self.depths['WellB5'][mini]}
                maxi = {"abscisse":0, "depth":self.depths['WellB5'][maxi]}

                # lines equations
                d = {}

                # for d1
                a = (float(m_coord["depth"]) - float(mini["depth"]))/(750)
                b = float(mini["depth"])
                d["d1"] = [a, b]

                # for d2
                a = alpha2/45
                d["d2"] = [a, b]

                # for d3
                a = - alpha1/45
                b = float(maxi["depth"])
                d["d3"] = [a, b]

                # for d4
                a = (float(m_coord["depth"]) - float(maxi["depth"]))/(750)
                d["d4"] = [a, b]

                if d["d3"][0] > d["d4"][0]:
                    print("Please choose a more important alpha angle")
                    return 0

                # Lines intersections

                Pts = {}

                # intersection d2-d3 : Point A
                xa = (d["d3"][1] - d["d2"][1])/(d["d2"][0] - d["d3"][0])
                ya = d["d2"][0] * xa + d["d2"][1]
                Pts["A"] = [xa, ya]

                # intersection d1-d3 : Point B
                xb = (d["d3"][1] - d["d1"][1])/(d["d1"][0] - d["d3"][0])
                yb = d["d1"][0] * xb + d["d1"][1]
                Pts["B"] = [xb, yb]

                # Point C : marker m
                xc = 750
                yc = float(m_coord["depth"])
                Pts["C"] = [xc, yc]

                # intersection d2-d4 : Point D
                xd = (d["d4"][1] - d["d2"][1])/(d["d2"][0] - d["d4"][0])
                yd = d["d2"][0] * xd + d["d2"][1]
                Pts["D"] = [xd, yd]

                ### discontinuity position

                ## Abscisse of the discontinuity
                # parameters of the normal law
                x1, x2 = max(xa, xb, xc, xd), min(xa, xb, xc, xd)
                x_mu = x2 + (x1 - x2)/2
                x_sigma = (x1 - x2)/6

                ## Depth of the discontinuity
                y1, y2 = max(ya, yb, yc, yd), min(ya, yb, yc, yd)
                y_mu = y2 + (y1 - y2) / 2
                y_sigma = (y1 - y2) / 6

                limits = [x1, x2, y1, y2]
                x_discontinuity, y_discontinuity = self.normal_law(x_mu, x_sigma, y_mu, y_sigma, n, limits, d)

                X, Y = np.array([xa, xb, xd, xc]), np.array([ya, yb, yd, yc])
                plt.plot(X, Y, 'x', color='b')
                plt.plot(x_discontinuity, y_discontinuity, 'o', color='r')
                plt.grid()
                plt.show()


    def biseau(self, well_name, marker_id, n, alpha1=None, alpha2=None):
        """
        :param alpha1:
        :param alpha2:
        :return: Position of all the discontinuities
        """
        if alpha1 == None:
            alpha1 = 22
        if alpha2 == None:
            alpha2 = 22
        w1 = self.wells_List[0]
        w2 = self.wells_List[1]

        if well_name == w1.well_name:
            m = w1.get_marker(marker_id)
            m_coord = {"abscisse": 0, "depth": self.depths['WellB5'][m.id]}
            correl = []
            for correl_m in m.correlations:
                correl.append(correl_m.id)
            mini, maxi = min(correl), max(correl)
            mini = {"abscisse": 750, "depth": self.depths['WellB11'][mini]}
            maxi = {"abscisse": 750, "depth": self.depths['WellB11'][maxi]}

            # lines equations
            d = {}

            # for d1
            print(m_coord["depth"])
            print(mini["depth"])
            print(maxi["depth"])
            a = (float(m_coord["depth"]) - float(mini["depth"])) / (-750)
            b = float(m_coord["depth"])
            d["d1"] = [a, b]

            # for d2
            a = - alpha2 / 45
            print(mini["depth"])
            b = float(mini["depth"]) - a * 750
            d["d2"] = [a, b]

            # for d3
            a = alpha1 / 45
            print(maxi["depth"])
            b = float(maxi["depth"]) - a * 750
            d["d3"] = [a, b]

            # for d4
            a = (float(m_coord["depth"]) - float(maxi["depth"])) / (-750)
            b = float(m_coord["depth"])
            d["d4"] = [a, b]
            print(d)

            if d["d3"][0] < d["d4"][0]:
                print("Please choose a more important alpha angle")
                return 0

            # Lines intersections

            Pts = {}

            # intersection d2-d3 : Point A
            xa = (d["d3"][1] - d["d2"][1]) / (d["d2"][0] - d["d3"][0])
            ya = d["d2"][0] * xa + d["d2"][1]
            Pts["A"] = [xa, ya]

            # intersection d1-d3 : Point B
            xb = (d["d3"][1] - d["d1"][1]) / (d["d1"][0] - d["d3"][0])
            yb = d["d1"][0] * xb + d["d1"][1]
            Pts["B"] = [xb, yb]

            # Point C : marker m
            xc = 0
            yc = float(m_coord["depth"])
            Pts["C"] = [xc, yc]

            # intersection d2-d4 : Point D
            xd = (d["d4"][1] - d["d2"][1]) / (d["d2"][0] - d["d4"][0])
            yd = d["d2"][0] * xd + d["d2"][1]
            Pts["D"] = [xd, yd]

            ### discontinuity position

            ## Abscisse of the discontinuity
            # parameters of the normal law
            x1, x2 = max(xa, xb, xc, xd), min(xa, xb, xc, xd)
            x_mu = x2 + (x1 - x2) / 2
            x_sigma = (x1 - x2) / 6

            ## Depth of the discontinuity
            y1, y2 = max(ya, yb, yc, yd), min(ya, yb, yc, yd)
            y_mu = y2 + (y1 - y2) / 2
            y_sigma = (y1 - y2) / 6

            limits = [x1, x2, y1, y2]
            x_discontinuity, y_discontinuity = self.normal_law(x_mu, x_sigma, y_mu, y_sigma, n, limits, d)

            X, Y = np.array([xa, xb, xd, xc]), np.array([ya, yb, yd, yc])
            plt.plot(X, Y, 'x', color='b')
            plt.plot(x_discontinuity, y_discontinuity, 'o', color='r')
            plt.grid()
            plt.show()



        else:
            m = w2.get_marker(marker_id)
            m_coord = {"abscisse": 750, "depth": self.depths['WellB11'][m.id]}
            correl = []
            for correl_m in m.correlations:
                correl.append(correl_m.id)
            mini, maxi = min(correl), max(correl)
            mini = {"abscisse": 0, "depth": self.depths['WellB5'][mini]}
            maxi = {"abscisse": 0, "depth": self.depths['WellB5'][maxi]}

            # lines equations
            d = {}

            # for d1
            a = (float(m_coord["depth"]) - float(mini["depth"])) / (750)
            b = float(mini["depth"])
            d["d1"] = [a, b]

            # for d2
            a = alpha2 / 45
            d["d2"] = [a, b]

            # for d3
            a = - alpha1 / 45
            b = float(maxi["depth"])
            d["d3"] = [a, b]

            # for d4
            a = (float(m_coord["depth"]) - float(maxi["depth"])) / (750)
            d["d4"] = [a, b]

            if d["d3"][0] > d["d4"][0]:
                print("Please choose a more important alpha angle")
                return 0

            # Lines intersections

            Pts = {}

            # intersection d2-d3 : Point A
            xa = (d["d3"][1] - d["d2"][1]) / (d["d2"][0] - d["d3"][0])
            ya = d["d2"][0] * xa + d["d2"][1]
            Pts["A"] = [xa, ya]

            # intersection d1-d3 : Point B
            xb = (d["d3"][1] - d["d1"][1]) / (d["d1"][0] - d["d3"][0])
            yb = d["d1"][0] * xb + d["d1"][1]
            Pts["B"] = [xb, yb]

            # Point C : marker m
            xc = 750
            yc = float(m_coord["depth"])
            Pts["C"] = [xc, yc]

            # intersection d2-d4 : Point D
            xd = (d["d4"][1] - d["d2"][1]) / (d["d2"][0] - d["d4"][0])
            yd = d["d2"][0] * xd + d["d2"][1]
            Pts["D"] = [xd, yd]


            ### discontinuity position

            ## Abscisse of the discontinuity
            # parameters of the normal law
            x1, x2 = max(xa, xb, xc, xd), min(xa, xb, xc, xd)
            x_mu = x2 + (x1 - x2) / 2
            x_sigma = (x1 - x2) / 6

            ## Depth of the discontinuity
            y1, y2 = max(ya, yb, yc, yd), min(ya, yb, yc, yd)
            y_mu = y2 + (y1 - y2) / 2
            y_sigma = (y1 - y2) / 6

            limits = [x1, x2, y1, y2]
            x_discontinuity, y_discontinuity = self.normal_law(x_mu, x_sigma, y_mu, y_sigma, n, limits, d)

            X, Y = np.array([xa, xb, xd, xc]), np.array([ya, yb, yd, yc])
            plt.plot(X, Y, 'x', color='b')
            plt.plot(x_discontinuity, y_discontinuity, 'o', color='r')
            plt.grid()
            plt.show()

    def normal_law(self, x_mu, x_sigma, y_mu, y_sigma, n, limits, d):
        X, Y = [], []
        xmax, xmin, ymax, ymin = limits

        for i in range(n):
            x = np.random.randn(1) * x_sigma + x_mu
            y = np.random.randn(1) * y_sigma + y_mu
            a1, b1, a2, b2, a3, b3, a4, b4 = d["d1"][0], d["d1"][1], d["d2"][0], d["d2"][1], d["d3"][0], d["d3"][1], \
                                             d["d4"][0], d["d4"][1]
            y1, y2, y3, y4 = x * a1 + b1, x * a2 + b2, x * a3 + b3, x * a4 + b4
            while not ((xmin < x[0] < xmax) and (ymin < y[0] < ymax) and (y > y1 and y > y3 and y < y2 and y < y4)):
                x = np.random.randn(1) * x_sigma + x_mu
                y = np.random.randn(1) * y_sigma + y_mu
                y1, y2, y3, y4 = x * a1 + b1, x * a2 + b2, x * a3 + b3, x * a4 + b4
            X.append(x[0])
            Y.append(y[0])
        print(X)
        print(Y)
        return X, Y


if __name__ == '__main__':
    # B = Biseau('result.csv')
    B = Biseau('testb5b11_wells.txt', 'result.csv').biseaux(100)
    # B = Biseau('test1_wells.txt', 'result_test.csv')
    # B = Biseau('testb5b11_wells.txt', 'result.csv').biseau("WellB5", 2, 100)
