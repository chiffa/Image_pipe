import csv
import os



def iterate_and_check_img(filename):
    f1 = open(filename, 'rb')
    f2 = open(filename+'.tmp', 'wb')
    writer = csv.writer(f2, delimiter = '\t')
    reader = csv.reader(f1, delimiter = '\t')

    for row in reader:
        if row[1] == '0':
            print row[0]
            row[1] = '1'
            writer.writerow(row)
        else:
            writer.writerow(row)

    f1.close()
    f2.close()
    os.rename(filename+'.tmp', filename)

    return filename

if __name__ == "__main__":
    iterate_and_check_img('test.tsv')
