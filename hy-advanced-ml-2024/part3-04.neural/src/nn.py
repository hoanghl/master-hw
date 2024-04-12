import sys
import numpy as np


class Edge:
    def __init__(self, inn, outn):
        """ 
        Creates an edge between two neurons (inn -> outn) with 0 weight

        Parameters
        ----------
        inn : Neuron
            Neuron sending its output
        outn : Neuron
            Neuron receveing the output
        """

        self.inn = inn
        self.outn = outn
        self.w = 0


class Neuron:
    def __init__(self):
        """ 
        Creates a neuron with no incoming or outgoing edges
        """

        self.inedges = []  # list of pointers for edge instances (=all incoming edges)
        self.bias = 0
        self.outedges = [] # list of pointers for edge instances (=all outgoing edges)

        # note that the same edge instance, say e, appears in e.inn.outedges _and_ in e.outn.inedges
        
        self.input = 0     # cached value before the activation function = weighted outputs of the inneurons + bias
        self.output = 0    # cached output = act(input)
        self.delta = 0     # cached delta value for doing backwards step


    def act(self, x):
        """ 
        Computes activation function (sigmoid) of x

        Parameters
        ----------
        x : real
            input

        Returns
        -------
        y : real
            sigmoid(x)
        """

        # place your code here
        return 1.0 / (1 + np.exp(-x))


    def der(self, x):
        """ 
        Computes the derivative of the activation function (sigmoid) of x

        Parameters
        ----------
        x : real
            input

        Returns
        -------
        y : real
            d sigmoid(x) / dx
        """

        # place your code here
        return self.act(x) * (1 - self.act(x))


    def compute_output(self):
        """ 
        Sets self.input and self.output based on the outputs of incoming neurons.
        self.input is the weighted sum of the incoming outputs + bias,
        self.output is the value after the activation function.
        """

        # compute self.input and self.output
        # place your code here
        self.input = 0
        for edge in self.inedges:
            self.input += edge.inn.output * edge.w 
        self.input += self.bias
        self.output = self.act(self.input)


    def compute_delta(self):
        """ 
        Computes self.delta based on the deltas of outgoing neurons and the derivative.
        """

        # compute self.delta
        # place your code here
        self.delta = 0
        for edge in self.outedges:
            self.delta += edge.w * edge.outn.delta
        self.delta *= self.der(self.input)


class NN:
    def __init__(self, ncnt):
        """ 
        Initializes a neural net with ncnt neurons and no edges.
    
        Parameters
        ----------
        ncnt: int
            number of neurons
        """

        self.neurons = [Neuron() for i in range(ncnt)] # list of neurons
        self.edges = [] # list of edges between neurons

    def join(self, i, j):
        """ 
        Adds an edge from ith to jth neuron. 

        Parameters
        ----------
        i : int
            index of a neuron
        j : int
            index of a neuron, must have i < j
        """

        assert(i < j)  # guarantees feed forward nn
        n1 = self.neurons[i]
        n2 = self.neurons[j]
        edge = Edge(n1, n2)
        n1.outedges.append(edge)
        n2.inedges.append(edge)
        self.edges.append(edge)


    def randomize(self):
        """ 
        Randomizes weights and edges (naive approach with no scaling tricks)
        """

        # randomize edge weights
        for e in self.edges:
            e.w = np.random.randn()

        # randomize biases
        for n in self.neurons:
            n.bias = np.random.randn()



    def forward(self, x):
        """ 
        Sets the outputs ith neuron to be x[i] for each 0 <= i < n, and
        computes the outputs of the remaining neurons using compute_output().
        Used for prediction and for forward phase in BP.

        Parameters
        ----------
        x : array of size n
            input values for the neural network
        """

        cnt = len(x)

        # set the output for the first layer
        for i in range(cnt):
            self.neurons[i].output = x[i]

        # compute the output for the remaining neurons
        # place your code here
        for i in range(cnt, len(self.neurons)):
            self.neurons[i].compute_output()


    def backward(self, y):
        """ 
        Computes the delta of the last neuron based on y. Then computes the
        deltas for the remaining neurons. Used in backward phase in BP.


        Parameters
        ----------
        y : real
            output of the training data, must be either 0 or 1
        """

        n = self.neurons[-1]
        n.delta = -n.der(n.input)
        if y != 1:
            n.delta = -n.delta

        # place your code here
        for i in range(len(self.neurons) - 2, -1, -1):
            self.neurons[i].compute_delta()
            

    def gradients(self, X, y):
        """ 
        Computes the gradients for the edge weights and for the bias terms as
        well as the error with the current weights.  The gradients and the
        error are normalized by the number of data points.

        Parameters
        ----------
        X : an array of size (n, k)
            n data points
        y : an array of size n
            n labels, values must be 0 or 1
    
        Returns
        -------
        wgrad : an array of ecnt
            wgrad[i] is the gradient for the weight of ith edge
        bgrad : an array of ncnt
            bgrad[i] is the gradient for the bias of ith neuron
        error : real
            Average error based on the output of the last neuron.  Note that
            the error is not the misclassication rate.
        """

        cnt = X.shape[0]       # number of data points

        ncnt = len(self.neurons)
        ecnt = len(self.edges)
        wgrad = np.zeros(ecnt) # gradient for weights
        bgrad = np.zeros(ncnt) # gradient for biases
        error = 0

        for i in range(cnt):
            self.forward(X[i, :])  # compute outputs
            self.backward(y[i])    # compute deltas

            # Update error
            error += np.abs(self.neurons[-1].output - y[i]) / cnt

            # Update the gradients 
            # place your code here
            for j in range(ncnt - 1, -1, -1):
                bgrad[j] += self.neurons[j].delta
            for j in range(ecnt - 1, -1, -1):
                wgrad[j] += self.edges[j].outn.delta * self.edges[j].inn.output

        bgrad /= cnt
        wgrad /= cnt

        return wgrad, bgrad, error


    def fit(self, X, y, itercnt, eta):
        """ 
        A naive implementation of BP algorithm.

        Parameters
        ----------
        X : an array of size (n, k)
            n data points
        y : an array of size n
            n labels, values must be 0 or 1
        itercnt : int
            number of iterations
        eta : real
            step size for the gradient
        """

        for i in range(itercnt):
            # compute the gradient
            wgrad, bgrad, error = self.gradients(X, y)
            print(error)

            # apply the gradient to edge weights
            for j, e in enumerate(self.edges):
                e.w -= eta*wgrad[j]

            # apply the gradient to biases
            for j, n in enumerate(self.neurons):
                n.bias -= eta*bgrad[j]


def main(argv):
    D = np.loadtxt(argv[1])

    incnt = D.shape[1] - 1
    hiddencnt = int(argv[2])
    itercnt = int(argv[3])
    eta = float(argv[4])

    y = D[:,0]
    X = D[:,1:]

    net = NN(incnt + hiddencnt + 1)

    for i in range(incnt):
        for j in range(hiddencnt):
            net.join(i, incnt + j)

    for i in range(hiddencnt):
        net.join(incnt + i, incnt + hiddencnt)
    

    # initial values for weights that do well for toy.txt
    #net.edges[1*hiddencnt + 0].w = 10
    #net.edges[2*hiddencnt + 1].w = 10

    #net.edges[incnt*hiddencnt + 0].w = -10
    #net.edges[incnt*hiddencnt + 1].w = -10
    #net.neurons[-1].bias = 1
    

    np.random.seed(100)
    net.randomize()

    net.fit(X, y, itercnt, eta)



# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__": 
    if len(sys.argv) != 5:
        print('usage: python %s filename hiddencnt itercnt eta' % sys.argv[0])
    else:
        main(sys.argv)
