#!/usr/bin/env python3

### This is a Python script I wrote during a DTU course on algorithms in bioinformatics
### It is an attempt at writing a Naïve Bayes classifier from scratch
### (as opposed to simply importing pre-written libraries)
### Code written by Jonas Dalsberg Jørgensen in June 2022

# Load modules
import numpy as np
import pandas as pd
from math import sqrt, pi, exp
from random import seed, randrange
import matplotlib.pyplot as plt
import itertools

# Load the cell type data format data
dataset = pd.read_csv("../../data/Project/cell_type_df.csv")
meta = dataset.pop("cell_type")
dataset.insert(dataset.shape[1], "cell_type", meta)
dataset = dataset.to_numpy().tolist()
class_cell = ["Th1", "Th17"]

for i in range(len(dataset)):
    if dataset[i][-1] == "th17":
        dataset[i][-1] = 1
    else:
        dataset[i][-1] = 0

# Load the treatment/control data (uncomment to load it)
treat_data = pd.read_csv("../../data/Project/combined_df.csv")
treat_data = treat_data.to_numpy()[:,1:].tolist()
class_treat = ["control", "treated"]

for i in range(len(treat_data)):
    if treat_data[i][-1] == "treated":
        treat_data[i][-1] = 1
    else:
        treat_data[i][-1] = 0

# Define functions
def separate_by_class(dataset):
    """" Split the dataset by class values, returns a dictionary """
    separated = dict()

    for i in range(len(dataset)):
        vector = dataset[i]
        class_value = vector[-1]

        if (class_value not in separated):
            separated[class_value] = list()

        separated[class_value].append(vector)

    return separated


def mean(numbers):
    """ Calculate the mean of a list of numbers """

    return sum(numbers)/float(len(numbers))


def stdev(numbers):
    """ Calculate the standard deviation of a list of numbers """
    avg = mean(numbers)
    variance = sum([(x-avg)**2 for x in numbers]) / float(len(numbers)-1)

    return sqrt(variance)


def summarize_dataset(dataset):
    """ Calculate the mean, stdev and count for each column in a dataset """
    summaries = [(mean(column), stdev(column), len(column)) for column in zip(*dataset)]
    del(summaries[-1])

    return summaries


def summarize_by_class(dataset):
    """ Split a dataset by class then calculate statistics for each row """
    separated = separate_by_class(dataset)
    summaries = dict()

    for class_value, rows in separated.items():
        summaries[class_value] = summarize_dataset(rows)

    return summaries


def calculate_probability(x, mean, stdev):
    """ Calculate the Gaussian probability distribution function for x """
    exponent = exp(-((x-mean)**2 / (2 * stdev**2 )))

    return (1 / (sqrt(2 * pi) * stdev)) * exponent


def calculate_class_probabilities(summaries, row):
    """ Calculate the probabilities of predicting each class for a given row """
    """ P(class=0|X1,X2) = P(X1|class=0) * P(X2|class=0) * P(class=0) """
    total_rows = sum([summaries[label][0][2] for label in summaries])
    probabilities = dict()

    for class_value, class_summaries in summaries.items():
        # probability of belonging to class "class_value", e.g. P(class = class_value)
        # (simply the no. of data points belonging to class e.g. 0 / the total no. of data points)
        probabilities[class_value] = summaries[class_value][0][2]/float(total_rows)

        # calculate the probability of having each gene expression value, e.g.
        # P(X1 | class = class_value), P(X2 | class = class_value), and so on
        # where X1,X2,X3 are gene expression values from a test observation
        for i in range(len(class_summaries)):
            mean, stdev, _ = class_summaries[i]
            probabilities[class_value] *= calculate_probability(row[i], mean, stdev)

    return probabilities


def predict(summaries, row):
    """ Predict the class for a given row """
    probabilities = calculate_class_probabilities(summaries, row)

    # initialise variables
    best_label, best_prob = None, -1

    for class_value, probability in probabilities.items():

        if best_label is None or probability > best_prob:
            best_prob = probability
            best_label = class_value

    return best_label


def naive_bayes(train, test):
    """ Naïve Bayes Algorithm """
    summaries = summarize_by_class(train)
    predictions = list()

    # for each test observation, predict the class
    for row in test:
        output = predict(summaries, row)
        predictions.append(output)

    return(predictions)


def cross_validation_split(dataset, n_folds):
    """ Split a dataset into k folds """
    dataset_split = list()
    dataset_copy = list(dataset)
    fold_size = int(len(dataset) / n_folds)

    for i in range(n_folds):
        fold = list()

        # append random data points until len(fold) == fold_size
        while len(fold) < fold_size:
            index = randrange(len(dataset_copy))
            fold.append(dataset_copy.pop(index))

        dataset_split.append(fold)

    return dataset_split


def accuracy_metric(actual, predicted):
    """ Calculate accuracy percentage """
    correct = 0

    for i in range(len(actual)):
        if actual[i] == predicted[i]:
            correct += 1

    return correct / float(len(actual)) * 100.0

def performance_measurement(actual, predicted):
    """ Calculates various performance measurements """
    TP = 0
    FP = 0
    TN = 0
    FN = 0

    for i in range(len(actual)):
        if actual[i] == predicted[i] == 1:
            TP += 1
        elif actual[i] == 0 and predicted[i] == 1:
            FP += 1
        elif actual[i] == predicted[i] == 0:
            TN += 1
        elif actual[i] == 1 and predicted[i] == 0:
            FN += 1

    # Sensitivity / recall / hit rate / true positive rate
    TPR = TP / (TP + FN) * 100

    # Specificity / true negative rate
    TNR = TN / (TN + FP) * 100

    # Precision / positive predictive value
    PPV = TP / (TP + FP) * 100

    return TP, FP, TN, FN, TPR, TNR, PPV

def evaluate_algorithm_avg(dataset, algorithm, n_folds):
    folds = cross_validation_split(dataset, n_folds)

    scores = list()
    TPs = list(); FPs = list(); TNs = list(); FNs = list()
    TPRs = list(); TNRs = list(); PPVs = list()

    for fold in folds:
        train_set = list(folds) # add all partitions to the training set
        train_set.remove(fold) # remove the test set from the training set
        train_set = sum(train_set, []) # flatten the list of lists
        test_set = list()

        # why not simply append the fold and flatten it?
        for row in fold:
            row_copy = list(row)
            test_set.append(row_copy)
            row_copy[-1] = None

        predicted = algorithm(train_set, test_set)
        actual = [row[-1] for row in fold]
        accuracy = accuracy_metric(actual, predicted)
        scores.append(accuracy)
        TP, FP, TN, FN, TPR, TNR, PPV = performance_measurement(actual, predicted)
        TPs.append(TP); FPs.append(FP); TNs.append(TN); FNs.append(FN)
        TPRs.append(TPR); TNRs.append(TNR); PPVs.append(PPV)

    return TPs, FPs, TNs, FNs, TPRs, TNRs, PPVs, scores

def evaluate_algorithm_cat(dataset, algorithm, n_folds):
    folds = cross_validation_split(dataset, n_folds)

    predicted = list()
    actual = list()

    for fold in folds:
        train_set = list(folds) # add all partitions to the training set
        train_set.remove(fold) # remove the test set from the training set
        train_set = sum(train_set, []) # flatten the list of lists
        test_set = list()

        # why not simply append the fold and flatten it?
        for row in fold:
            row_copy = list(row)
            test_set.append(row_copy)
            row_copy[-1] = None

        predicted += algorithm(train_set, test_set)
        actual += [row[-1] for row in fold]

    TP, FP, TN, FN, TPR, TNR, PPV = performance_measurement(actual, predicted)
    cnf_matrix = np.array(([[TN, FP],[FN,TP]]), dtype=np.int64)
    accuracy = accuracy_metric(actual, predicted)

    return cnf_matrix, TPR, TNR, PPV, accuracy

def plot_confusion_matrix(cnf_matrix, class_names, title="Confusion matrix"):
    """ Plots a confusion matrix """
    plt.figure(figsize = (9,6))
    plt.imshow(cnf_matrix, interpolation="nearest", cmap=plt.get_cmap("Blues"))
    plt.title(title)
    plt.colorbar()

    tick_marks = np.arange(len(class_names))
    plt.xticks(tick_marks, class_names, rotation=45)
    plt.yticks(tick_marks, class_names)

    fmt = "d"
    thresh = cnf_matrix.max() / 2.

    for i, j in itertools.product(range(cnf_matrix.shape[0]), range(cnf_matrix.shape[1])):
        plt.text(j, i, format(cnf_matrix[i, j], fmt), horizontalalignment = "center", color="white" if cnf_matrix[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel("True label")
    plt.xlabel("Predicted label")

def multi_cnf_matrices(TNs, FPs, FNs, TPs):
    """ Produces multiple confusion matrices from lists of TPs, FPs, FNs, and TNs """
    cnf_matrices = list()
    for i in range(len(TNs)):
        cnf_matrices.append(np.array(([[TNs[i], FPs[i]],[FNs[i], TPs[i]]]), dtype=np.int64))

    return cnf_matrices

def subplot_cnf(cnf_matrix, class_names, plot_no, no_plots = 5):
    """ Makes a confusion matrix designed to be used by the function plot_cnf_multiple """
    u = int(np.floor(np.sqrt(no_plots))); v = int(np.ceil(float(no_plots)/u)) # determine how many plots to have per row and column
    plt.subplot(u, v, plot_no)
    plt.imshow(cnf_matrix, interpolation="nearest", cmap=plt.get_cmap("Blues"))
    #alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    plt.colorbar()

    tick_marks = np.arange(len(class_names))
    plt.xticks(tick_marks, class_names, rotation=45)
    plt.yticks(tick_marks, class_names)

    fmt = "d"
    thresh = cnf_matrix.max() / 2.

    for i, j in itertools.product(range(cnf_matrix.shape[0]), range(cnf_matrix.shape[1])):
        plt.text(j, i, format(cnf_matrix[i, j], fmt), horizontalalignment = "center", color="white" if cnf_matrix[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel("True label")
    plt.xlabel("Predicted label")
    plt.title("Fold {}".format(plot_no), loc="left", y=1.1)

def plot_cnf_multiple(multi_cnfs, class_names, title="Confusion matrix", no_plots = 5):
    """ Plots a confusion matrix """
    plt.figure(figsize = (14,7))
    plt.suptitle(title, fontsize = 20, y=1)
    for i in range(no_plots):
        subplot_cnf(multi_cnfs[i], class_names, plot_no = i+1, no_plots = no_plots)


# Cross-validation loop
rand_seed = 1
print("Naïve Bayes classifier performance:")
print("Th1 vs Th17")
seed(rand_seed)
TPs, FPs, TNs, FNs, TPRs, TNRs, PPVs, acc_cell = evaluate_algorithm_avg(dataset, naive_bayes, 5)
cnf_multi_cell = multi_cnf_matrices(TNs, FPs, FNs, TPs)
print("TP: {}\tFP: {}\tTN: {}\tFN: {}".format(sum(TPs), sum(FPs), sum(TNs), sum(FNs)))

print("\nPerformance per fold:")
print("Sensitivity / true positive rate:\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%".format(TPRs[0], TPRs[1], TPRs[2], TPRs[3], TPRs[4]))
print("Specificity / true negative rate:\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%".format(TNRs[0], TNRs[1], TNRs[2], TNRs[3], TNRs[4]))
print("Precision:\t\t\t\t\t\t\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%".format(PPVs[0], PPVs[1], PPVs[2], PPVs[3], PPVs[4]))
print("Accuracy:\t\t\t\t\t\t\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%".format(acc_cell[0], acc_cell[1], acc_cell[2], acc_cell[3], acc_cell[4]))

# print("\nAveraging over the folds:")
# print("Sensitivity / true positive rate:\t{:.2f}%".format(sum(TPRs)/float(len(TPRs))))
# print("Specificity / true negative rate:\t{:.2f}%".format(sum(TNRs)/float(len(TNRs))))
# print("Precision:\t\t\t\t\t\t\t{:.2f}%".format(sum(PPVs)/float(len(PPVs))))
# print("Accuracy:\t\t\t\t\t\t\t{:.2f}%".format(sum(acc_cell)/float(len(acc_cell))))

print("\nConcatenating the folds:")
seed(rand_seed)
cnf_cell, TPR, TNR, PPV, accuracy = evaluate_algorithm_cat(dataset, naive_bayes, 5)
print("Sensitivity / true positive rate:\t{:.2f}%".format(TPR))
print("Specificity / true negative rate:\t{:.2f}%".format(TNR))
print("Precision:\t\t\t\t\t\t\t{:.2f}%".format(PPV))
print("Accuracy:\t\t\t\t\t\t\t{:.2f}%".format(accuracy))


print("\n\nControl vs treatment")
seed(rand_seed)
TPs, FPs, TNs, FNs, TPRs, TNRs, PPVs, acc_treat = evaluate_algorithm_avg(treat_data, naive_bayes, 5)
cnf_multi_treat = multi_cnf_matrices(TNs, FPs, FNs, TPs)
print("TP: {}\tFP: {}\tTN: {}\tFN: {}".format(sum(TPs), sum(FPs), sum(TNs), sum(FNs)))

print("\nPerformance per fold:")
print("Sensitivity / true positive rate:\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%".format(TPRs[0], TPRs[1], TPRs[2], TPRs[3], TPRs[4]))
print("Specificity / true negative rate:\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%".format(TNRs[0], TNRs[1], TNRs[2], TNRs[3], TNRs[4]))
print("Precision:\t\t\t\t\t\t\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%".format(PPVs[0], PPVs[1], PPVs[2], PPVs[3], PPVs[4]))
print("Accuracy:\t\t\t\t\t\t\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%\t{:.2f}%".format(acc_treat[0], acc_treat[1], acc_treat[2], acc_treat[3], acc_treat[4]))


# # print("\nAveraging over the folds:")
# # print("Sensitivity / true positive rate:\t{:.2f}%".format(sum(TPRs)/float(len(TPRs))))
# # print("Specificity / true negative rate:\t{:.2f}%".format(sum(TNRs)/float(len(TNRs))))
# # print("Precision:\t\t\t\t\t\t\t{:.2f}%".format(sum(PPVs)/float(len(PPVs))))
# # print("Accuracy:\t\t\t\t\t\t\t{:.2f}%".format(sum(acc_treat)/float(len(acc_treat))))

print("\nConcatenating the folds:")
seed(rand_seed)
cnf_treat, TPR, TNR, PPV, accuracy = evaluate_algorithm_cat(treat_data, naive_bayes, 5)
print("Sensitivity / true positive rate:\t{:.2f}%".format(TPR))
print("Specificity / true negative rate:\t{:.2f}%".format(TNR))
print("Precision:\t\t\t\t\t\t\t{:.2f}%".format(PPV))
print("Accuracy:\t\t\t\t\t\t\t{:.2f}%".format(accuracy))


# Plotting the accuracies of the two classifiers
X = np.arange(5)
fig,ax = plt.subplots(figsize=(8,3))
test_cell = ax.bar(X, acc_cell, width = 0.35, label =  'Th1 vs Th17')
test_treat = ax.bar(X + 0.35, acc_treat, width = 0.35, label =  'Control vs treatment')
plt.ylabel("Accuracy", fontsize=10);
ax.set_xticks([0.175, 1.175, 2.175, 3.175, 4.175])
ax.set_xticklabels(["1", "2", "3", "4", "5"])
plt.ylim(0, 100)
plt.xlabel("Models", fontsize=10);
plt.title("Accuracy of NB classifier models")
plt.legend(loc = 'lower right')
plt.tight_layout()
plt.savefig("../../plots/NB_accuracy_barplot.png", dpi = 100, bbox_inches = "tight")

# #Plotting the confusion matrix
plot_confusion_matrix(cnf_cell, class_cell)
plt.savefig("../../plots/cnf_cell.png", dpi = 100, bbox_inches = "tight")
# plot_confusion_matrix(cnf_treat, class_treat)


# # Multiple confusion matrices plotted
# plot_cnf_multiple(cnf_multi_treat, class_treat, title = "Confusion matrix for each fold (control vs. treatment)")
# plt.savefig("../../plots/cnf_multi_treat.png", dpi = 100, bbox_inches = "tight")

# plot_cnf_multiple(cnf_multi_cell, class_cell, title = "Confusion matrix for each fold (Th1 vs. Th17)")
# plt.savefig("../../plots/cnf_multi_cell.png", dpi = 100, bbox_inches = "tight")
