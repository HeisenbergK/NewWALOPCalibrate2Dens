import numpy as np


# This is a function that takes: ([[q,u],[q,u],...], [a0,a0,...], [aq,aq,...], [au,au,...], [aq2,aq2,...], [au2,au2,
# ...], [aqu,aqu,...]) and returns the equivalent polynomial
def modell(inpolcoords, a0, aq, au, aq2, au2, aqu):
    newcoords = list(map(list, zip(*inpolcoords)))
    qs = newcoords[0]
    us = newcoords[1]
    return np.multiply(a0, np.ones_like(qs)) + \
        np.multiply(aq, qs) + np.multiply(au, us) + \
        np.multiply(aq2, np.multiply(qs, qs)) + \
        np.multiply(au2, np.multiply(us, us)) + \
        np.multiply(aqu, np.multiply(qs, us))


# This is a function that takes: ([q,u], a0, aq, au, aq2, au2, aqu) and returns the equivalent polynomial
def modellrun(inpolcoords, a0, aq, au, aq2, au2, aqu):
    qs = inpolcoords[0]
    us = inpolcoords[1]
    return a0 + (aq * qs) + (au * us) + (aq2 * qs * qs) + (au2 * us * us) + (aqu * qs * us)


def modellnm(inpolcoords, a0, aq, au, aq2, au2):
    newcoords = list(map(list, zip(*inpolcoords)))
    qs = newcoords[0]
    us = newcoords[1]
    return np.multiply(a0, np.ones_like(qs)) + \
        np.multiply(aq, qs) + np.multiply(au, us) + \
        np.multiply(aq2, np.multiply(qs, qs)) + \
        np.multiply(au2, np.multiply(us, us))


# This is a function that takes: ([q,u], a0, aq, au, aq2, au2, aqu) and returns the equivalent polynomial
def modellnmrun(inpolcoords, a0, aq, au, aq2, au2):
    qs = inpolcoords[0]
    us = inpolcoords[1]
    return a0 + (aq * qs) + (au * us) + (aq2 * qs * qs) + (au2 * us * us)


def modell3(inpolcoords, a0, aq, au, aq2, au2, aqu, aq3, au3, aq2u, aqu2):
    newcoords = list(map(list, zip(*inpolcoords)))
    qs = newcoords[0]
    us = newcoords[1]
    return np.multiply(a0, np.ones_like(qs)) + \
        np.multiply(aq, qs) + np.multiply(au, us) + \
        np.multiply(aq2, np.multiply(qs, qs)) + \
        np.multiply(au2, np.multiply(us, us)) + \
        np.multiply(aqu, np.multiply(qs, us)) + \
        np.multiply(aq3, np.multiply(qs, np.multiply(qs, qs))) + \
        np.multiply(au3, np.multiply(us, np.multiply(us, us))) + \
        np.multiply(aq2u, np.multiply(qs, np.multiply(qs, us))) + \
        np.multiply(aqu2, np.multiply(qs, np.multiply(us, us)))


# This is a function that takes: ([q,u], a0, aq, au, aq2, au2, aqu) and returns the equivalent polynomial
def modell3run(inpolcoords, a0, aq, au, aq2, au2, aqu, aq3, au3, aq2u, aqu2):
    qs = inpolcoords[0]
    us = inpolcoords[1]
    return a0 + \
        (aq * qs) + \
        (au * us) + \
        (aq2 * qs * qs) + \
        (au2 * us * us) + \
        (aqu * qs * us) + \
        (aq3 * qs * qs * qs) + \
        (au3 * us * us * us) + \
        (aq2u * qs * qs * us) + \
        (aqu2 * qs * us * us)

# This is a function that takes: ([[q,u],[q,u],...], [a0,a0,...], [aq,aq,...], [au,au,...], [aq2,aq2,...], [au2,au2,
# ...], [aqu,aqu,...]) and returns the equivalent polynomial


def modellconst(inpolcoords, aq, au, aq2, au2, aqu):
    newcoords = list(map(list, zip(*inpolcoords)))
    qs = newcoords[0]
    us = newcoords[1]
    a0 = newcoords[2]
    return np.multiply(a0, np.ones_like(qs)) + \
        np.multiply(aq, qs) + np.multiply(au, us) + \
        np.multiply(aq2, np.multiply(qs, qs)) + \
        np.multiply(au2, np.multiply(us, us)) + \
        np.multiply(aqu, np.multiply(qs, us))


# This is a function that takes: ([q,u], a0, aq, au, aq2, au2, aqu) and returns the equivalent polynomial
def modellrunconst(inpolcoords, aq, au, aq2, au2, aqu):
    qs = inpolcoords[0]
    us = inpolcoords[1]
    a0 = inpolcoords[2]
    return a0 + (aq * qs) + (au * us) + (aq2 * qs * qs) + (au2 * us * us) + (aqu * qs * us)


def modellnmconst(inpolcoords, aq, au, aq2, au2):
    newcoords = list(map(list, zip(*inpolcoords)))
    qs = newcoords[0]
    us = newcoords[1]
    a0 = newcoords[2]
    return np.multiply(a0, np.ones_like(qs)) + \
        np.multiply(aq, qs) + np.multiply(au, us) + \
        np.multiply(aq2, np.multiply(qs, qs)) + \
        np.multiply(au2, np.multiply(us, us))


# This is a function that takes: ([q,u], a0, aq, au, aq2, au2, aqu) and returns the equivalent polynomial
def modellnmrunconst(inpolcoords, aq, au, aq2, au2):
    qs = inpolcoords[0]
    us = inpolcoords[1]
    a0 = inpolcoords[2]
    return a0 + (aq * qs) + (au * us) + (aq2 * qs * qs) + (au2 * us * us)


def modell3const(inpolcoords, aq, au, aq2, au2, aqu, aq3, au3, aq2u, aqu2):
    newcoords = list(map(list, zip(*inpolcoords)))
    qs = newcoords[0]
    us = newcoords[1]
    a0 = newcoords[2]
    return np.multiply(a0, np.ones_like(qs)) + \
        np.multiply(aq, qs) + np.multiply(au, us) + \
        np.multiply(aq2, np.multiply(qs, qs)) + \
        np.multiply(au2, np.multiply(us, us)) + \
        np.multiply(aqu, np.multiply(qs, us)) + \
        np.multiply(aq3, np.multiply(qs, np.multiply(qs, qs))) + \
        np.multiply(au3, np.multiply(us, np.multiply(us, us))) + \
        np.multiply(aq2u, np.multiply(qs, np.multiply(qs, us))) + \
        np.multiply(aqu2, np.multiply(qs, np.multiply(us, us)))


# This is a function that takes: ([q,u], a0, aq, au, aq2, au2, aqu) and returns the equivalent polynomial
def modell3runconst(inpolcoords, aq, au, aq2, au2, aqu, aq3, au3, aq2u, aqu2):
    qs = inpolcoords[0]
    us = inpolcoords[1]
    a0 = inpolcoords[2]
    return a0 + \
        (aq * qs) + \
        (au * us) + \
        (aq2 * qs * qs) + \
        (au2 * us * us) + \
        (aqu * qs * us) + \
        (aq3 * qs * qs * qs) + \
        (au3 * us * us * us) + \
        (aq2u * qs * qs * us) + \
        (aqu2 * qs * us * us)
