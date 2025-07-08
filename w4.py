import numpy as np
import json
from matplotlib import pyplot as plt


class w4_project:
    def __init__(self, name):

        self.name = name  # the input json file name
        self.alpha_L0 = 0  # thin airfoil zero lift angle of attack = 0

        self.airfoil_lift_slope = None
        self.nodes_per_semispan = None
        self.n = None

        # planform parameters
        self.type = None
        self.aspect_ratio = None
        self.taper_ratio = None
        self.filename = None

        # washout parameters
        self.distribution = None
        self.amount_deg = None
        self.CL_design = None
        self.B3 = None

        # aileron parameters
        self.begin_zb = None
        self.end_zb = None
        self.begin_cfc = None
        self.end_cfc = None
        self.hinge_efficiency = None
        self.deflection_efficiency = 1

        # wing conditions
        self.alpha_root_deg = None
        self.CL = None
        self.aileron_deflection_deg = None
        self.pbar = None

        # view
        self.view_planform = None
        self.view_washout = None
        self.view_aileron = None
        self.view_CL_hat = None
        self.view_CL_tilde = None

        # point data from file
        self.c_b_source = None
        self.z_b_source = None

        # calculated variables
        self.theta = None
        self.c_b = None
        self.z_b = None
        self.C_matrix = None
        self.C_inv = None
        self.an = None
        self.bn = None
        self.omega = None  # distribution
        self.CL_alpha = None
        self.Omega = None  # washout amount
        self.eomega = None
        self.alpha_root = None
        self.CL1 = None
        self.CL2 = None
        self.CDi1 = None
        self.CDi2 = None
        self.cf_c = None
        self.chi = None
        self.cn = None
        self.dn = None
        self.Clpbar = None
        self.Clda = None
        self.An = None
        self.CDI = None
        self.Cl1 = None
        self.Cl2 = None

        # CL hat plot parameters
        self.CL_hat = None
        self.CL_hat_planform = None
        self.CL_hat_washout = None
        self.CL_hat_aileron = None
        self.CL_hat_roll = None

        # CL tide plot parameters
        self.CL_tilde = None
        self.CL_tilde_planform = None
        self.CL_tilde_washout = None
        self.CL_tilde_aileron = None
        self.CL_tilde_roll = None

    def open(self):
        # load parameters from json file

        json_string = open(self.name).read()
        input_dict = json.loads(json_string)

        self.airfoil_lift_slope = input_dict["wing"]["airfoil_lift_slope"]
        self.nodes_per_semispan = input_dict["wing"]["nodes_per_semispan"]
        self.n = int(2*self.nodes_per_semispan-1)

        # planform
        self.type = input_dict["wing"]["planform"]["type"]
        self.aspect_ratio = input_dict["wing"]["planform"]["aspect_ratio"]
        self.taper_ratio = input_dict["wing"]["planform"]["taper_ratio"]
        self.filename = input_dict["wing"]["planform"]["filename"]

        # washout
        self.distribution = input_dict["wing"]["washout"]["distribution"]
        self.amount_deg = input_dict["wing"]["washout"]["amount[deg]"]
        self.CL_design = input_dict["wing"]["washout"]["CL_design"]
        self.B3 = input_dict["wing"]["washout"]["B3"]

        # aileron
        self.begin_zb = input_dict["wing"]["aileron"]["begin[z/b]"]
        self.end_zb = input_dict["wing"]["aileron"]["end[z/b]"]
        self.begin_cfc = input_dict["wing"]["aileron"]["begin[cf/c]"]
        self.end_cfc = input_dict["wing"]["aileron"]["end[cf/c]"]
        self.hinge_efficiency = input_dict["wing"]["aileron"]["hinge_efficiency"]

        # condition
        self.alpha_root_deg = input_dict["condition"]["alpha_root[deg]"]
        self.CL = input_dict["condition"]["CL"]
        self.aileron_deflection_deg = input_dict["condition"]["aileron_deflection[deg]"]
        self.pbar = input_dict["condition"]["pbar"]

        # view
        self.view_planform = input_dict["view"]["planform"]
        self.view_washout = input_dict["view"]["washout_distribution"]
        self.view_aileron = input_dict["view"]["aileron_distribution"]
        self.view_CL_hat = input_dict["view"]["CL_hat_distributions"]
        self.view_CL_tilde = input_dict["view"]["CL_tilde_distributions"]

    def get_theta(self):
        # get theta
        self.theta = np.linspace(0.0, np.pi, self.n)

    def get_cb_zb(self):
        # get z/b
        self.z_b = -0.5*np.cos(self.theta)

        # get c/b when planform is a file
        if self.type == "file":
            # pull out the point and remove c/b and z/b text
            with open(self.filename, "r") as f:
                points = []
                c_b = []
                z_b = []
                for line in f:
                    points = line.split()
                    c_b.append(points[1])
                    z_b.append(points[0])
            f.close()
            c_b.remove("c/b")
            z_b.remove("z/b")

            # recalculcate c/b with z/b from theta
            self.c_b_source = [float(coordinate) for coordinate in c_b]
            self.z_b_source = [float(coordinate) for coordinate in z_b]
            self.c_b = np.interp(
                abs(self.z_b), self.z_b_source, self.c_b_source)
            wingsurface = 0
            for i in range(0, len(c_b)-1):
                wingsurface = (self.c_b_source[i] + self.c_b_source[i+1])*abs(
                    self.z_b_source[i]-self.z_b_source[i+1])*0.5 + wingsurface
            self.aspect_ratio = (2*self.z_b_source[-1])**2/(2*wingsurface)

        # get c/b when planform is a elliptic wing or tapered wing
        else:
            if self.type == "elliptic":
                # (6.27)
                self.c_b = 4*np.sin(self.theta)/(np.pi*self.aspect_ratio)
            if self.type == "tapered":
                # (6.28)
                self.c_b = 2/(self.aspect_ratio*(1+self.taper_ratio)) * \
                    (1-(1-self.taper_ratio)*np.abs(np.cos(self.theta)))

    def get_C_matrix(self):
        # calculate C matrix
        n = self.n
        zeroindex = np.argwhere(self.c_b == 0)

        # if there is a c/b at wingtip replace it with 0.001
        for i in zeroindex:
            self.c_b[i] = 10**(-3)

        # calculate C matrix
        self.C_matrix = np.zeros(shape=(n, n))
        for i in range(0, n):
            for j in range(1, n+1):
                if i == 0:
                    self.C_matrix[i][j-1] = (j)**2
                elif i == n-1:
                    self.C_matrix[i][j-1] = (-1)**(j+1)*(j)**2
                else:
                    self.C_matrix[i][j-1] = (4/(self.c_b[i]*self.airfoil_lift_slope) + (
                        j) / np.sin(self.theta[i]))*(np.sin((j)*self.theta[i]))

    def get_C_inv(self):
        # get C inv
        id_matrix = np.eye(self.n)
        self.C_inv = np.linalg.solve(self.C_matrix, id_matrix)

    def get_distribution(self):
        # calculcate washout distribution (omega)
        theta = self.theta
        n = self.n
        zeroindex = np.argwhere(theta == 0)
        for i in zeroindex:
            theta[i] = 10**(-3)

        if self.distribution == "none":
            self.omega = np.zeros(shape=(n, 1))

        if self.distribution == "linear":
            # eqution (6.30)
            self.omega = np.zeros(shape=(n, 1))
            for i in range(0, n):
                self.omega[i] = abs(np.cos(theta[i]))

        if self.distribution == "optimum":
            # eqation (6.38)
            self.omega = np.zeros(shape=(n, 1))
            for i in range(0, n):
                P1 = (4/self.airfoil_lift_slope)
                P2 = (1-self.B3)/self.c_b[self.nodes_per_semispan-1]
                P3 = (np.sin(self.theta[i]) + self.B3 *
                      np.sin(3*self.theta[i]))/self.c_b[i]
                P4 = 3*self.B3 * \
                    (1+np.sin(3*self.theta[i]) / np.sin(self.theta[i]))
                P5 = 4*(1-self.B3)/(self.airfoil_lift_slope *
                                    self.c_b[self.nodes_per_semispan-1]) - 12*self.B3
                self.omega[i] = (P1*(P2-P3)-P4)/P5

    def get_an_bn(self):
        ones = np.ones(shape=(self.n, 1))
        # get an
        self.an = np.dot(self.C_inv, ones)
        # get bn
        self.bn = np.dot(self.C_inv, self.omega)

    def get_K(self):
        # get kappa parameters

        # get kappa_L (6.21)
        self.KL = (1-(1+np.pi*self.aspect_ratio/self.airfoil_lift_slope) *
                   self.an[0])/((1+np.pi*self.aspect_ratio/self.airfoil_lift_slope)*self.an[0])

        # get kappa_D (6.16)
        self.KD = 0
        for j in range(2, self.n+1):
            self.KD = self.KD + (self.an[j-1]/self.an[0])**2*j

        # get kappa_DL (6.17)
        self.KDL = 0
        for j in range(2, self.n+1):
            self.KDL = self.KDL + 2*self.bn[0]/self.an[0] * \
                (self.an[j-1]/self.an[0])*j * \
                (self.bn[j-1]/self.bn[0] - self.an[j-1]/self.an[0])
        if self.distribution == "none":
            self.KDL = np.array([0])

        # get kappa_D_Omega (6.18)
        self.KDomega = 0
        for j in range(2, self.n+1):
            self.KDomega = self.KDomega + ((self.bn[0]/self.an[0])**2) * \
                j*((self.bn[j-1]/self.bn[0] - self.an[j-1]/self.an[0])**2)
        if self.distribution == "none":
            self.KDomega = np.array([0])

    def get_CL_alpha(self):
        # get CL_alpha (6.20)
        self.CL_alpha = self.airfoil_lift_slope / \
            ((1+self.airfoil_lift_slope/(np.pi*self.aspect_ratio))*(1+self.KL))

    def get_washoutamount(self):
        # get washout amount (Omega)

        # if there is a number for washout amount in json file
        if self.amount_deg != "optimum":
            self.Omega = np.deg2rad(self.amount_deg)

        # if the washout is optimum
        if self.amount_deg == "optimum":

            if self.distribution == "optimum":
                # equation (6.39)
                self.Omega = self.CL_design/np.pi/self.aspect_ratio * \
                    (4*(1-self.B3)/self.airfoil_lift_slope /
                     self.c_b[self.nodes_per_semispan-1] - 12*self.B3)

            elif self.type == "elliptic":
                # washout amount = 0 for elliptic wing
                self.Omega = np.array([0])

            else:
                # use (6.32) for the rest case
                self.Omega = self.KDL*self.CL_design / \
                    (2*self.KDomega*self.CL_alpha)

        if self.distribution == "none":
            # no washout Omega = 0
            self.Omega = np.array([0])

    def get_eomega(self):
        # get epilson_omega with (6.22)
        self.eomega = self.bn[0]/self.an[0]

    def get_alpha_root(self):
        # use (6.19) if use CL
        if self.alpha_root_deg == "CL":
            self.alpha_root = self.CL/self.CL_alpha + \
                self.eomega*self.Omega + self.alpha_L0
        else:
            # use value from the json file
            self.alpha_root = np.deg2rad(self.alpha_root_deg)

    def get_cfc(self):
        # find aileron x coordinate
        aileron_x = np.array(
            [-self.end_zb, -self.begin_zb, self.begin_zb, self.end_zb])

        # find aileron x coordinate
        if self.type == "elliptic":
            # find y using (6.27)
            aileron_y = 4/(np.pi*self.aspect_ratio) * \
                (1-(2*aileron_x)**2)**0.5  # (6.27)

        elif self.type == "tapered":
            # find y using (6.28)
            aileron_y = 2 / \
                (self.aspect_ratio*(1+self.taper_ratio)) * \
                (1-(1-self.taper_ratio)*abs(2*aileron_x))  # (6.28)
        else:
            # find y if planform is a file
            aileron_y = np.interp(
                abs(aileron_x), self.z_b_source, self.c_b_source)

        # just a array make calculation easiler
        c = np.array([self.end_cfc, self.begin_cfc,
                      self.begin_cfc, self.end_cfc])

        # compute the slope of aileron leading edge
        m1 = ((-0.75+c[0]) * aileron_y[0]-(-0.75+c[1])
              * aileron_y[1])/(aileron_x[0]-aileron_x[1])
        m2 = ((-0.75+c[2]) * aileron_y[2]-(-0.75+c[3]) *
              aileron_y[3])/(aileron_x[2] - aileron_x[3])

        # compute cf_c using slope
        cf_c = np.zeros(shape=(self.n, 1))
        for i in range(0, self.n):
            if self.begin_zb <= self.z_b[i] <= self.end_zb:
                cf_c[i] = ((-0.75+c[2])*aileron_y[2] + m2 *
                           (self.z_b[i]-self.begin_zb) - (-0.75*self.c_b[i]))/self.c_b[i]
            if -self.end_zb <= self.z_b[i] <= -self.begin_zb:
                cf_c[i] = ((-0.75+c[0])*aileron_y[0] + m1 *
                           (self.z_b[i]-(-self.end_zb)) - (-0.75*self.c_b[i]))/self.c_b[i]
        self.cf_c = cf_c

    def get_chi(self):
        # get aileron distribution
        # (6.36)
        costheta = 2*self.cf_c-1
        theta_f = np.arccos(2*self.cf_c-1)

        # (6.35)
        epsilon_fi = 1 - (theta_f - (1 - costheta**2)**0.5)/np.pi

        chi = np.zeros(shape=(self.n, 1))
        # epsilon_f (6.34)
        for i in range(0, self.n):
            if i > self.n*0.5:
                # right side aileron distribution
                chi[i] = -self.hinge_efficiency * \
                    self.deflection_efficiency*epsilon_fi[i]
            else:
                # left side aileron distribution
                chi[i] = self.hinge_efficiency * \
                    self.deflection_efficiency*epsilon_fi[i]

        self.chi = chi

    def get_cn_dn(self):
        # get cn & dn array
        self.cn = np.dot(self.C_inv, self.chi)
        self.dn = np.dot(self.C_inv, np.cos(self.theta)).reshape(self.n, 1)
        # reshape just to make cn and dn at the same form(2D n*1 array)

    def get_pbar(self):
        # get pbar
        if self.pbar != "steady":
            self.Clpbar = 0
            self.Clda = 0
        if self.pbar == "steady":
            # (6.24)
            self.Clda = - np.pi*self.aspect_ratio*self.cn[1]/4
            # (6.25)
            self.Clpbar = - np.pi*self.aspect_ratio*self.dn[1]/4
            # (6.26)
            self.pbar = -self.Clda/self.Clpbar * \
                np.deg2rad(self.aileron_deflection_deg)

    def get_An(self):
        # get An using (6.10)
        An = np.zeros(shape=(self.n, 1))
        for i in range(0, self.n):
            An[i] = self.an[i]*(self.alpha_root-self.alpha_L0) - self.bn[i] * \
                self.Omega + \
                self.cn[i]*np.deg2rad(self.aileron_deflection_deg) + \
                self.dn[i]*self.pbar
        self.An = An

    def get_CL(self):

        # form (6.5)
        self.CL1 = np.pi*self.aspect_ratio*self.An[0]

        # from (6.19)
        self.CL2 = self.CL_alpha * \
            ((self.alpha_root - self.alpha_L0) - self.eomega*self.Omega)

    def get_CDi(self):
        # (6.6) includes effects of ailerons and rolling rate
        CDI = 0
        for j in range(0, self.n):
            CDI = CDI + np.pi*self.aspect_ratio * (j+1)*self.An[j]**2
        self.CDi1 = CDI - np.pi*self.aspect_ratio*self.pbar*self.An[1]/2

        # (6.15) neglects effects of ailerons and rolling rate
        self.CDi2 = (self.CL2**2*(1+self.KD) - self.KDL*self.CL2*self.CL_alpha*self.Omega + self.KDomega *
                     (self.CL_alpha*self.Omega)**2)/(np.pi*self.aspect_ratio)

    def get_Cl(self):
        # get rolling moment
        # (6.7)
        self.Cl1 = -np.pi * self.aspect_ratio*self.An[1]/4

        # (6.23)
        self.Cl2 = self.Clda * \
            np.deg2rad(self.aileron_deflection_deg) + \
            self.Clpbar*self.pbar

    def get_Cn(self):
        Cn = 0

        # get Cn with (6.9)
        for j in range(1, self.n):
            Cn = Cn + np.pi*self.aspect_ratio/4 * \
                (2*(j+1)-1)*self.An[j-1]*self.An[j]
        self.Cn = Cn - np.pi*self.aspect_ratio * \
            self.pbar/8*(self.An[0]+self.An[2])

    def show(self):
        # display the result
        print("------- Solution --------------------------")
        print("aspect ratio = ", self.aspect_ratio)
        print("------- Kappa Results ---------------------")
        print("KL = ", self.KL)
        print("CL_alpha = ", self.CL_alpha)
        print("eomega = ", self.eomega)
        print("KD = ", self.KD)
        print("KDL = ", self.KDL)
        print("KDomega = ", self.KDomega)
        print("Clda", self.Clda)
        print("Clpbar", self.Clpbar)
        print("------- Operating Condition Results -------")
        print("alpha_root_deg = ", np.rad2deg(self.alpha_root))
        # print("alpha_root = ", self.alpha_root) # in radius
        print("Omega_deg = ", np.rad2deg(self.Omega))
        # print("Omega = ", self.Omega) # in radius
        print("aileron = ", self.aileron_deflection_deg)
        print("(6.5) CL = ", self.CL1)
        print("(6.19) CL = ", self.CL2)
        print("(6.6) CDi = ", self.CDi1)
        print("(6.15) CDi = ", self.CDi2)
        # print("(6.7) Cl = ", self.Cl1) # Cl for equation 6.7
        print("(6.23) Cl = ", self.Cl2)
        print("Cn = ", self.Cn)
        print("pbar = ", self.pbar)
        print("-------------------------------------------")

    def plot_wing(self):
        if self.view_planform == True:
            x = self.z_b[:]
            y = 0.25*self.c_b[:]
            y1 = -0.75*self.c_b[:]
            ax = plt.subplot()

            # Plot leading edge & trailing edge & 1/4 chord
            ax.plot(x, y, color="black")
            ax.plot(x, y1, color="black")
            ax.plot(x, 0*y, color="black")
            plt.xlabel("z/b")
            plt.ylabel("c/b")
            plt.title("planform")

            # plot the blue section line
            for i in range(0, self.n):
                x2 = np.array([self.z_b[i], self.z_b[i]])
                y2 = np.array([0.25*self.c_b[i], -0.75 *
                              self.c_b[i] + self.c_b[i]*self.cf_c[i, 0]])
                ax.plot(x2, y2, color="blue")

            # get aileron parameter
            aileron_x = np.array(
                [-self.end_zb, -self.begin_zb, self.begin_zb, self.end_zb])
            if self.type == "elliptic":
                aileron_y = 4/(np.pi*self.aspect_ratio) * \
                    (1-(2*aileron_x)**2)**0.5
            elif self.type == "tapered":
                aileron_y = 2 / \
                    (self.aspect_ratio*(1+self.taper_ratio)) * \
                    (1-(1-self.taper_ratio)*abs(2*aileron_x))
            else:
                aileron_y = np.interp(
                    abs(aileron_x), self.z_b_source, self.c_b_source)

            # plot the aileron (red)
            for i in range(0, self.n):
                x6 = np.array([self.z_b[i], self.z_b[i]])
                y6 = np.array([-0.75*self.c_b[i] + self.c_b[i]
                              * self.cf_c[i, 0], -0.75*self.c_b[i]])
                ax.plot(x6, y6, color="red")

            c = np.array([self.end_cfc, self.begin_cfc,
                          self.begin_cfc, self.end_cfc])

            # Plot the begin and end aileron (black)
            for i in range(0, len(c)):
                x3 = np.array([aileron_x[i], aileron_x[i]])
                y3 = np.array([-0.75*aileron_y[i], (-0.75+c[i])*aileron_y[i]])
                ax.plot(x3, y3, color="black")
            x4 = np.array([aileron_x[0], aileron_x[1]])
            y4 = np.array([(-0.75+c[0])*aileron_y[0],
                          (-0.75+c[1])*aileron_y[1]])
            ax.plot(x4, y4, color="black")

            x5 = np.array([aileron_x[2], aileron_x[3]])
            y5 = np.array([(-0.75+c[2])*aileron_y[2],
                          (-0.75+c[3])*aileron_y[3]])
            ax.plot(x5, y5, color="black")

            ax.set_aspect("equal", adjustable="box")
            plt.show()

    def plot_distribution(self):
        # plot the washout distribution
        if self.view_washout == True:
            x = self.z_b[:]
            y = self.omega[:]
            bx = plt.subplot()
            bx.plot(x, y, color="black")
            plt.xlabel("z/b")
            plt.ylabel("omega")
            plt.title("Washout Distribution")
            plt.show()

    def plot_ailerondistribution(self):
        # plot the aileron distribution
        if self.view_aileron == True:
            x = self.z_b[:]
            y = self.chi[:]
            cx = plt.subplot()
            cx.plot(x, y, color="black")
            plt.xlabel("z/b")
            plt.ylabel("chi")
            plt.title("Aileron Distribution")
            plt.show()

    def plot_CL_hat(self):
        # plot CL_hat
        if self.view_CL_hat == True:
            # get CL_hat_planform (6.45)
            self.CL_hat_planform = np.zeros(self.n)
            for i in range(0, self.n):
                for j in range(0, self.n):
                    self.CL_hat_planform[i] += 4*(self.alpha_root-self.alpha_L0) * \
                        self.an[j]*np.sin((j+1)*self.theta[i])

            # get CL_hat_washout (6.46)
            self.CL_hat_washout = np.zeros(self.n)
            for i in range(0, self.n):
                for j in range(0, self.n):
                    self.CL_hat_washout[i] += -4*self.Omega * \
                        self.bn[j]*np.sin((j+1)*self.theta[i])

            # get CL_hat_aileron (6.47)
            self.CL_hat_aileron = np.zeros(self.n)
            for i in range(0, self.n):
                for j in range(0, self.n):
                    self.CL_hat_aileron[i] += 4*np.deg2rad(
                        self.aileron_deflection_deg)*self.cn[j]*np.sin((j+1)*self.theta[i])

            # get CL_hat_roll (6.48)
            self.CL_hat_roll = np.zeros(self.n)
            for i in range(0, self.n):
                for j in range(0, self.n):
                    self.CL_hat_roll[i] += 4*self.pbar * \
                        self.dn[j]*np.sin((j+1)*self.theta[i])

            # CL_hat (6.44)
            self.CL_hat = self.CL_hat_planform + self.CL_hat_washout + \
                self.CL_hat_aileron + self.CL_hat_roll

            # get coordinate
            x = self.z_b[:]
            y1 = self.CL_hat_planform[:]
            y2 = self.CL_hat_washout[:]
            y3 = self.CL_hat_aileron[:]
            y4 = self.CL_hat_roll[:]
            y5 = self.CL_hat[:]

            # plot
            dx = plt.subplot()
            dx.plot(x, y1, color="blue", label="planform")
            dx.plot(x, y2, color="red", label="washout")
            dx.plot(x, y3, color="green", label="aileron")
            dx.plot(x, y4, color="gray", label="roll rate")
            dx.plot(x, y5, color="black", label="total")
            dx.legend(loc="upper right")
            plt.xlabel("z/b")
            plt.ylabel("CL_hat")
            plt.title("CL_hat Distribution")
            plt.show()

    def plot_CL_tilde(self):
        if self.view_CL_tilde == True:
            # get CL_tilde_planform (6.51)
            self.CL_tilde_planform = np.zeros(self.n)
            for i in range(0, self.n):
                self.CL_tilde_planform[i] = self.CL_hat_planform[i]/self.c_b[i]

            # get CL_tilde_washout (6.52)
            self.CL_tilde_washout = np.zeros(self.n)
            for i in range(0, self.n):
                self.CL_tilde_washout[i] = self.CL_hat_washout[i]/self.c_b[i]

            # get CL_tilde_aileron (6.53)
            self.CL_tilde_aileron = np.zeros(self.n)
            for i in range(0, self.n):
                self.CL_tilde_aileron[i] = self.CL_hat_aileron[i]/self.c_b[i]

            # get CL_tilde_roll (6.54)
            self.CL_tilde_roll = np.zeros(self.n)
            for i in range(0, self.n):
                self.CL_tilde_roll[i] = self.CL_hat_roll[i]/self.c_b[i]

            # get CL_tilde (6.50)
            self.CL_tilde = self.CL_tilde_planform + self.CL_tilde_washout + \
                self.CL_tilde_aileron + self.CL_tilde_roll

            # get coordinate
            x = self.z_b[:]
            y1 = self.CL_tilde_planform[:]
            y2 = self.CL_tilde_washout[:]
            y3 = self.CL_tilde_aileron[:]
            y4 = self.CL_tilde_roll[:]
            y5 = self.CL_tilde[:]

            # plot
            ex = plt.subplot()
            ex.plot(x, y1, color="blue", label="planform")
            ex.plot(x, y2, color="red", label="washout")
            ex.plot(x, y3, color="green", label="aileron")
            ex.plot(x, y4, color="gray", label="roll rate")
            ex.plot(x, y5, color="black", label="total")
            ex.legend(loc="upper right")
            plt.xlabel("z/b")
            plt.ylabel("CL_hat")
            plt.title("CL_hat Distribution")

            plt.show()

    def get_bending_moment_distribution(self):
        # Reserved for extra points, but can't figure it out
        # unused
        C_mb = "somethings"

    def write_for_check(self):
        # To make it easier to check anwsers, the calculation results are written in separate TXT files
        # unused
        n = self.n
        filename = "aileron.txt"
        with open(filename, "w") as f:
            f.write("alieron\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.chi[i, 0]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "cfc.txt"
        with open(filename, "w") as f:
            f.write("cfc\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.cf_c[i, 0]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "An.txt"
        with open(filename, "w") as f:
            f.write("An\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.An[i, 0]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "an.txt"
        with open(filename, "w") as f:
            f.write("an\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.an[i, 0]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "bn.txt"
        with open(filename, "w") as f:
            f.write("bn\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.bn[i, 0]), end="\n", file=f)
            f.write("\n")

        filename = "cn.txt"
        with open(filename, "w") as f:
            f.write("cn\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.cn[i, 0]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "dn.txt"
        with open(filename, "w") as f:
            f.write("dn\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.dn[i, 0]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "omega.txt"
        with open(filename, "w") as f:
            f.write("omega\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.omega[i, 0]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "z_b.txt"
        with open(filename, "w") as f:
            f.write("z/b\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.z_b[i]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "c_b.txt"
        with open(filename, "w") as f:
            f.write("c/b\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.c_b[i]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "theta.txt"
        with open(filename, "w") as f:
            f.write("theta\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.theta[i]), end="\n", file=f)
            f.write("\n")
        f.close()

        filename = "sin.txt"
        with open(filename, "w") as f:
            f.write("sin\n")
            for i in range(0, n):
                print("{:< 17.12f}".format(
                    self.sin[i]), end="\n", file=f)
            f.write("\n")
        f.close()

    def write(self):
        # write the solutions.txt file
        n = self.n
        filename = "solutions.txt"
        with open(filename, "w") as f:

            # write C_matrix
            f.write("C_matrix\n")
            for i in range(0, n):
                for j in range(0, n):
                    print("{:< 17.12f}".format(
                        self.C_matrix[i, j]), end="", file=f)
                print("", end="\n", file=f)
            f.write("\n")

            # write C_inv
            f.write("C_inv\n")
            for i in range(0, n):
                for j in range(0, n):
                    print("{:< 17.12f}".format(
                        self.C_inv[i, j]), end="", file=f)
                print("", end="\n", file=f)
            f.write("\n")

            # write an
            f.write("an\n")
            for i in range(0, n):
                print("{:> 17.12f}".format(self.an[i, 0]), end="\n", file=f)
            f.write("\n")

            # write bn
            f.write("bn\n")
            for i in range(0, n):
                print("{:> 17.12f}".format(self.bn[i, 0]), end="\n", file=f)
            f.write("\n")

            f.write("cn\n")
            for i in range(0, n):
                print("{:> 17.12f}".format(self.cn[i, 0]), end="\n", file=f)
            f.write("\n")

            f.write("dn\n")
            for i in range(0, n):
                print("{:> 17.12f}".format(self.dn[i, 0]), end="\n", file=f)
            f.write("\n")

        f.close()

    def run(self):
        # run all functions needed
        self.open()
        self.get_theta()
        self.get_cb_zb()
        self.get_C_matrix()
        self.get_C_inv()
        self.get_distribution()
        self.get_an_bn()
        self.get_K()
        self.get_CL_alpha()
        self.get_washoutamount()
        self.get_eomega()
        self.get_alpha_root()
        self.get_cfc()
        self.get_chi()
        self.get_cn_dn()
        self.get_pbar()
        self.get_An()
        self.get_CL()
        self.get_CDi()
        self.get_Cl()
        self.get_Cn()
        self.show()
        self.plot_wing()
        self.plot_distribution()
        self.plot_ailerondistribution()
        self.plot_CL_hat()
        self.plot_CL_tilde()
        self.write()


if __name__ == "__main__":
    w4 = w4_project("input2.json")
    w4.run()
