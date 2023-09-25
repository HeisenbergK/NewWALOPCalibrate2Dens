import numpy as np


def qu(fieldx, fieldy):
    theta = np.arctan2(fieldy, fieldx)
    q = np.cos(np.multiply(2.0, theta))
    u = np.sin(np.multiply(2.0, theta))
    return q, u


def EVPAdeg(q,u):
    return (180.0/np.pi)*EVPArad(q,u)


def EVPArad(q,u):
    return (0.5*np.arctan2(u,q))+(np.pi/2.0)


def pol2cart(phi):
    x = np.cos(phi)
    y = np.sin(phi)
    return x, y


def trans_to_pol(t1,t2,t3,t4):
    if not (np.size(t1)==np.size(t2) and np.size(t2)==np.size(t3) and  np.size(t2)==np.size(t4)):
        raise Exception("Uneven arguments in trans_to_pol")
    q = np.divide(np.subtract(t1,t3),np.add(t1,t3))
    u = np.divide(np.subtract(t4,t2),np.add(t2,t4))
    return [q,u]