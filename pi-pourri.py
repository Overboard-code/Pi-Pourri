#!/usr/bin/env python3
# Pi-pourri.py
# By Andrew Richter
# March 2022
#
#   I started with picrunch.py by Don Cross
# See his wonderful article at:
# https://medium.com/@cosinekitty/how-to-calculate-a-million-digits-of-pi-d62ce3db8f58
#   I had a kid ask me "How long would it take to calcuate a Million digits of Pi?"
#   This code can Use Machin's Formula or one of the Machin Like formulae with gmpy2 mpz() for the numbers
#    Machin 1706
#  π/4 = 4*(4*arctan(1/5) - arctan(1/239)) or
#    Hwang Chien-Lih, 2003
#  π/4 = 183*arctan(1/239) + 32*arctan(1/1023)-68*arctan(5832)+
#         12*arctan(1/110443)-12*arctan(4841182)-100*arctan(6826318)
#    Hwang Chien-Lih, 2003
#  π/4 =  183*arctan(1/239) + 32*arctan(1/1023)-68*arctan(1/5832)+
#         12*arctan(1/113021)-100*arctan(1/6826318)-12*arctan(1/33366019650)+
#         12*arctan(1/43599522992503626068)
#  or Chudnovosky 1988: π = (Q(0, N) / 12T(0, N) + 12AQ(0, N))**(C**(3/2)
#
#   Try --help to see the list of available formulae. All have been checked using -d 1,000,000 against https://www.piday.org/million/
#   I also used the last five digits of known results to check the answers for lengths that are powers of ten (100,....1,000,000,000)
#   If you want to add a Machin-like_formula, just add it to setOfNames, setOfMults etc.  Each list item Must contain the same
#     number of entries (except name of course)
#   For Machin-like formulae each calulation for arctan(1/nnnnn) gets its own multiprocessing.Pool() thread.
#   When all the threads are done the rest of the formula is processed
#   It calculates pi to -d xxxx places after the decimal The Default is 100,000
#   and writes the answer to a -f filename The default is pi.txt in the current directory
#
#   All of the formulae can get to a million digits in a few seconds  Real differences start to show up around 10 million.
#   It takes about two and a half minutes for Chudnovsky(8) to get to 100 Million digits
#   The parser will allow -d as high as a Hundred Million. Chudnovsky and AGM fail at a Billion because th calulation is too large for GMPY.
#   100 million digits takes anywhere from 3 to 45 minutes depending on which formula is used
#   I used crude timeing using time() - start_time to generate elapsed seconds  There are better ways
#
import sys,time,multiprocessing,unicodedata,logging,os,argparse
from datetime import timedelta
from functools import partial
try:
    from gmpy2 import mpz,mpq,mul,isqrt,mpfr,atan2,tan,atan,sqrt,get_context  # Gumpy2 mpz large ints are ten times faster than python large int
except ImportError:
    raise ImportError('This program requires gmpy2. exiting....')
    sys.exit(1)
# Change logging to INFO or WARNING to see less output
logging.basicConfig(level=("DEBUG"),format='[%(levelname)s] %(asctime)s %(funcName)s: %(processName)s %(message)s')
#
def say_formula(credit,mults,denoms,signs):
    # format a Manchin like formula string from 3 lists and an author's credit
    saypi = unicodedata.lookup("GREEK SMALL LETTER PI")
    if ('Chudnovsky' in credit) or ('AGM' in credit) or ("Bellard" in credit):
        form = credit
    else:
        form = "\t{}\n\t ".format(credit) + saypi + "/4 = "
        for i in range(0,len(denoms)):
            if i == 0:
                sign = '' # No leading sign of first arctan multiple
            else:
                sign = '\t\t' + (('-','+')[signs[i]>0])  +  ' '
            if mults[i] == 1: # Don't say 1*arctan() just say arctan()
                mx = ""
            else:
                mx = str(mults[i]) + "*"
            form = form + ' ' + sign + mx + 'arctan(1/' +str(denoms[i]) + ')\n'
    return form

class PiAGM:
    LOG2_10 = 3.321928094887362

    def __init__(self,ndigits):
        self.ndigits = ndigits
        self.cdigits = ndigits + len(str(ndigits))+9      # Extra digits to reduce trailing error More factors means more error
        self.iters = 0

    def compute(self):
        get_context().precision=int(self.cdigits * self.LOG2_10)
        epsilon = mpfr(1/mpfr(10**self.ndigits))
        logging.debug('AGM precision({:,}) Started '.format(self.ndigits ) )
        start_time = time.time()
        a = mpfr(1)
        b = mpfr(1/sqrt(mpfr(2)))

        diff = mpfr(a - b)
        series = n = mpfr (0)
        while diff > epsilon:
            series += mpfr(2)**(n) * (diff**mpfr(2))
            n += mpfr(1)
            arith = mpfr((a + b)/mpfr (2))
            geom = mpfr(sqrt(a*b))
            a, b = arith, geom
            diff = mpfr(a - b)
            self.iters += 1
            if (self.iters % 10  == 0 ):
                logging.debug('AGM ... {:,} iterations and {:.2f} seconds.'.format( int(self.iters),time.time() - start_time))
        # a and b have converged to the AGM
        get_context().precision=int((self.ndigits+2 ) * self.LOG2_10)
        pi = mpfr(4)*a*a/(mpfr(1) - series)
        logging.debug('AGM Done! {:,} iterations and {:.2f} seconds.'.format(self.iters,time.time() - start_time) )
        return  str(pi)[:-2],self.iters,time.time()-start_time

class PiBellard:
        LOG2_10 = 3.321928094887362
        def __init__(self,ndigits):
            self.ndigits = ndigits
            self.iters = 0

        def compute(self):
            cdigits = self.ndigits+15
            get_context().precision=int(cdigits * self.LOG2_10) # Precision isn't digits  need some math
            start_time = time.time()  # Start the clock for total time
            logging.debug('Bellard precision({:,}) Started '.format(self.ndigits ) )
             #http://en.wikipedia.org/wiki/Bellard%27s_formula
            pi = mpfr(0)
            for i in range(self.ndigits):
                a = mpfr(1)/(16**i)
                b = mpfr(4)/(8*i+1)
                c = mpfr(2)/(8*i+4)
                d = mpfr(1)/(8*i+5)
                e = mpfr(1)/(8*i+6)
                r = mpfr(a*(b-c-d-e))
                pi += r
                self.iters += 1
                if (self.iters % 10000  == 0 ):
                    logging.debug('Bellard ... {:,} iterations and {:.2f} seconds.'.format( int(self.iters),time.time() - start_time))
            get_context().precision=int((self.ndigits+4) * self.LOG2_10) # Precision isn't digits  need some math
            pi = pi + 0
            logging.debug('Bellard Done! {:,} iterations and {:.2f} seconds.'.format(self.iters,time.time() - start_time) )
            return str(pi)[:self.ndigits+2],self.iters,time.time()-start_time


class PiMachin:
    LOG2_10 = 3.321928094887362
    def __init__(self,ndigits,name,denoms,mults,operators):
        """ Initialization
        :param int digits: digits of PI computation
        :param string name: name for credit on formula
        :param list  denoms a lis of ints for the denomiators
        :param list  mults: a list of ints as multipliyers for the Machin formula
        :param list operators: a list of 1 or -1 to cause addition or subtraction
        """
        self.name = name
        self.ndigits = ndigits
        self.denoms = denoms
        self.mults = mults
        self.operators = operators
        self.xdigits = len(denoms)+7              # Extra digits to reduce trailing error More factors means more error

    def ArctanDenom(self,d):

        cdigits = self.ndigits+self.xdigits
        get_context().precision=int(cdigits * self.LOG2_10)
        # Calculates arctan(1/d) = 1/d - 1/(3*d^3) + 1/(5*d^5) - 1/(7*d^7) + ...
        logging.debug('arctan(1/{} Started ) '.format(d ) )
        total = mpfr()
        arc_start_time = time.time()  # Start the clock for this arctan calulation
        total = mpfr(atan2(mpfr(1),mpfr(d)))
        logging.debug('arctan(1/{}) Done!   {:.2f} seconds.'.format(int(d),time.time() - arc_start_time))
        return total,int(1)
    #
    def compute(self):
        start_time = time.time()  # Start the clock for total time
        name = self.name  # pull the chosen formula list from the list of formulae
        denoms = self.denoms
        mults = self.mults
        operators = self.operators
        ndigits = self.ndigits
        cdigits = self.ndigits+self.xdigits
        logging.info("Starting Machin-Like formula to {:,} decimal places".format(ndigits) )
        get_context().precision=int(cdigits * self.LOG2_10)
        logging.debug("Starting {} Pool threads to calculate arctan values.".format(len(denoms)) )
        p =  multiprocessing.Pool(processes=(len(denoms))) # get some threads for our pool
        results=p.map(self.ArctanDenom, denoms) # one arctan(1/xxxx) per thread
        p.close()
        p.join()  # wait for them to finish
        # Now we have the arctan calculations from the pool threads in results[]
        # Apply chosen Formula to the results and calculate pi using mults and signs
        logging.debug ("Now multiplying and summing arctan results")
        arctanSum = pi = mpfr(0)
        i = iters = 0
        for tup in results:
            iters += tup[1]  # Keep track of this thread's iterations for later
            arctanSum += mpfr(mpfr(self.mults[i])*mpfr(tup[0])*mpfr(self.operators[i])) # Add or subtract the product from the accumulated arctans
            i += 1 # Next result
        get_context().precision=int((self.ndigits+2) * self.LOG2_10)
        pi = mpfr(4) * arctanSum # change pi/4 = x to pi = 4 * x
        # We calculated extra digits to compensate for roundoff error.
        # Chop off the extra digits now.
        return str(pi)[:self.ndigits+2],iters,time.time()-start_time

class PiChudnovsky:
    A = 13591409
    B = 545140134
    C = 640320
    D = 426880
    E = 10005
    C3_24  = C ** 3 // 24
    LOG2_10 = 3.321928094887362
    #DIGITS_PER_TERM = math.log(53360 ** 3) / math.log(10)  #=> 14.181647462725476
    DIGITS_PER_TERM = 14.181647462725476
    MMILL = mpz(1000000)


    def __init__(self,ndigits):
        """ Initialization
        :param int digits: digits of PI computation
        :param string filename: output file to write digits
        """
        self.ndigits = ndigits
        self.n      = mpz(self.ndigits // self.DIGITS_PER_TERM + 1)
        self.prec   = mpz((self.ndigits + 1) * self.LOG2_10)
        self.one_sq = mpz(10) ** (2 * ndigits)
        self.sqrt_c = isqrt(self.E * self.one_sq)
        self.iters  = 0

    def compute(self):
        """ Computation """
        try:
            start_time = time.time()
            logging.debug("Starting {} formula to {:,} decimal places".format(name,ndigits) )
            p, q, t = self.__bsa(0, self.n)
            pi = (q * self.D * self.sqrt_c) // t
            logging.debug('{} calulation Done! {:,} iterations and {:.2f} seconds.'.format( name, int(self.iters),time.time() - start_time))
            return str(pi),self.iters,time.time() - start_time
        except Exception as e:
            raise

    def __bsa(self, a, b):
        """ PQT computation by BSA(= Binary Splitting Algorithm)
        :param int a: positive integer
        :param int b: positive integer
        :return list [int p_ab, int q_ab, int t_ab]
        """
        try:
            self.iters += 1
            if (self.iters % self.MMILL  == mpz(0) ):
                logging.debug('Chudnovsky ... {:,} iterations and {:.2f} seconds.'.format( int(self.iters),time.time() - start_time))
            if a + 1 == b:
                if a == 0:
                    p_ab = q_ab = mpz(1)
                else:
                    p_ab = mpz((6 * a -5) * (2 * a - 1) * (6 * a - 1))
                    q_ab = mpz(a * a * a * self.C3_24)
                t_ab = p_ab * (self.A + self.B * a)
                if a & 1:
                    t_ab *= -1
            else:
                m = (a + b) // 2
                p_am, q_am, t_am = self.__bsa(a, m)
                p_mb, q_mb, t_mb = self.__bsa(m, b)
                p_ab = p_am * p_mb
                q_ab = q_am * q_mb
                t_ab = q_mb * t_am + p_am * t_mb
            return [p_ab, q_ab, t_ab]
        except Exception as e:
            raise


if __name__ == '__main__':
    last5DigitsOfPi = {
             10 : 26535,
            100 : 70679,
           1000 : 1989,
          10000 : 75678,
         100000 : 24646,
        1000000 : 58151,
        1234567 : 14707,
       10000000 : 55897,
      100000000 : 51592,
     1000000000 : 45519,
    }
#  Took values from lists from Machin and Miachin like formulae here:
#  https://en.wikipedia.org/wiki/Machin-like_formula
#  Add to all 4 lists to add a new formula  Each formula list entry *Must* be the same size for each type of list
    setOfNames = ["John Machin 1706",
        "F. C. M. Störmer 1896",
        "Kikuo Takano 1982",
        "Hwang Chien-Lih, 1997",
        "Hwang Chien-Lih, 2003",
        "Jörg Uwe Arndt 1993 ",
        "Hwang Chien-Lih, 2004",
        "\tRadius Generator- Fabrice Bellard?, 1997 \n\tπ = 126N∑n=0(−1)n210n(−254n+1−14n+3+2810n+1−2610n+3−2210n+5−2210n+7+110n+9)\n",
        "\tThe Square AGM - Salamin & Brent, 1976\n\tπ = limit as n goes to infinity  (an+bn)**2/(4tn)\n",
         "\tChudnovsky brothers  1988 \n\tπ = (Q(0, N) / 12T(0, N) + 12AQ(0, N))**(C**(3/2))\n" ]
    setOfMults = [ [4,1],  # Machin 1706
        [44,7,12,24], #  F. C. M. Störmer 1896
        [12,32,5,12], # Kikuo Takano (1982)
        [183,32,68,12,12,100], # Hwang 1997
        [183,32,68,12,100,12,12], # Hwang 2003
        [36462,135908,274509,39581,178477,114569,146571,61914,69044,89431,43938], # Jörg Uwe Arndt 1993
        [36462,26522,19275,3119,3833,5183,37185,11010,3880,16507,7476],  # Hwang 2004
        ["Place","holder"], # Improved Hex - Fabrice Bellard, 1997
        ["Place","holder"], # The Square AGM - Salamin & Brent, 1976
        ["Place","holder"] ]  # Chudnovsky brothers  1988
    setOfDenoms = [ [5,239],
        [57,239,682,12943],
        [49,57,239,110443],
        [239,1023,5832,110443,4841182,6826318],
        [239,1023,5832,113021,6826318,33366019650,43599522992503626068],
        [390112,485298,683982,1984933,2478328,3449051,18975991,22709274,24208144,201229582,2189376182],
        [51387,485298,683982,1984933,2478328,3449051,18975991,22709274,24208144,201229582,2189376182],
        ["Place","holder"], # Improved Hex - Fabrice Bellard, 1997
        ["Place","holder"], # The Square AGM - Salamin & Brent, 1976
        ["Place","holder"] ] # Chudnovsky brothers  1988
    setOfOpers = [ [1,-1],
        [1,1,-1,1],
        [1,1,-1,1],
        [1,1,-1,1,-1,-1],
        [1,1,-1,1,-1,-1,1],
        [1,1,1,-1,1,-1,-1,1,-1,-1,-1],
        [1,1,1,-1,-1,-1,-1,-1,1,-1,-1],
        ["Place","holder"], # Improved Hex/Radius - Fabrice Bellard, 1997
        ["Place","holder"], # The Square AGM - Salamin & Brent, 1976
        ["Place","holder"] ] # Chudnovsky brothers  1988

    numFormulae= len(setOfDenoms)
    formRange="[1 to {}]".format(numFormulae)
    def range_type(test_value, min=1, max=10):
        value = int(''.join(filter(lambda i: i.isdigit(), test_value)) )
        if min <= value <= max:
            return value
        else:
            raise argparse.ArgumentTypeError('value {} not in range {:,} to {:,}'.format(value,min,max))

    # create strings for each formula for help.  Add here and to the descStr if you added a formula to the lists above
    formList = ""
    for i in range(numFormulae):
        formList += " {} {} \n".format(i+1,say_formula(setOfNames[i],setOfMults[i],setOfDenoms[i],setOfOpers[i]))

    pgmName =  os.path.basename(sys.argv[0])
    descStr = """ {0} runs an algoritym from a list to calulate Pi to a number of decimal places
      Default: {0} --digits 100000 --file pi.txt --alog 4

      So -d 100,000,000 will take a while to finish, -d 1,000,000 very quickly
      A last 5 digit check is done on powers of ten (10,...10000000,100000000)
 eg.  {0} --file elbow.txt -d 1000000 -a 4
      {0} -f test.txt -d 123,456

      List of Formulae:

{1} """.format(pgmName,formList)
    parser = argparse.ArgumentParser(description=descStr,formatter_class=argparse.RawDescriptionHelpFormatter )
    # add expected arguments
    parser.add_argument('-f','--file', nargs='?', dest='filename', default='pi.txt',
                required=False,  help="File Name to write Pi to.. Default is %(default)s")
    parser.add_argument('-d','--digits', nargs=1, dest='max_digits', metavar="[1 to 100,000,000]", default=[100000],
                type=partial(range_type, min=1, max=100000000), required=False,
                help="How many digits to calculate.  Default is %(default)s ")
    parser.add_argument('-a','--algo',nargs=1, dest='algo', metavar=formRange, default=[4],
                type=partial(range_type, min=1, max=numFormulae), required=False, help="Which Machin(like) formula. Default is %(default)s")

    args = parser.parse_args(sys.argv[1:])
    if args.max_digits:
        ndigits = int(args.max_digits[0])
    if args.filename:
        outFileName =  args.filename
    if args.algo:
        algox =  args.algo[0] - 1

    start_time = time.time()  # Start the clock for total time

    name = setOfNames[algox]  # pull the chosen formula list from the list of formulae
    denoms = setOfDenoms[algox]
    mults = setOfMults[algox]
    operators = setOfOpers[algox]
    saypi = unicodedata.lookup("GREEK SMALL LETTER PI")

    logging.info("Computing {} to ( {:,} digits )".format(saypi,ndigits))

    if 'Chudnovsky' in name:
        obj = PiChudnovsky(ndigits)
    else:
        if 'AGM' in name:
            obj = PiAGM(ndigits)
        else:
            if 'Bellard' in name:
                obj = PiBellard(ndigits)
            else:
                obj = PiMachin(ndigits,name,denoms,mults,operators)

    pi,iters,tm = obj.compute()
    endDigits = int(pi[-5:])
    if ndigits in last5DigitsOfPi:
            if last5DigitsOfPi[ndigits] == endDigits:
                logging.info("Last 5 digits of {} were {} as expected at offset {:,}".format(saypi, endDigits,ndigits-5 ))
            else:
                logging.warning("\nWRONG WRONG WRONG")
                logging.warning("\nLast 5 digits were {} and are WRONG should be {}\n".format(endDigits, last5DigitsOfPi[ndigits]))
                logging.warning("\nWRONG WRONG WRONG")
    else:
        logging.info("Last five digits are {} but length {:,} wasn't in the list of known values".format(endDigits,ndigits))

    startWrite = time.time()
    with open(outFileName, 'wt') as outfile:
        outfile.write(pi)
    tw = time.time() - startWrite
    logging.info("Calculated {} to {:,} digits using a formula of:\n {} {} "
        .format(saypi,ndigits,algox+1,say_formula(name,mults,denoms,operators) ) )
    logging.info('Wrote {:,} digits of {} to file {} in {}'.format(ndigits,saypi,outFileName,str(timedelta(seconds=tw))))
    logging.info("Calculation and write took: {:,} iterations and  {}."
        .format(int(iters),str(timedelta(seconds=tm))) )
    sys.exit(0)
