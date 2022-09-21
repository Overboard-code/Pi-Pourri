#!/usr/bin/env python3
# Pi-pourri.py
# By Andrew Richter
# March 2022
#
#   I started with picrunch.py by Don Cross
# See his wonderful article at:
# https://medium.com/@cosinekitty/how-to-calculate-a-million-digits-of-pi-d62ce3db8f58
# and See: http://www.craig-wood.com/nick/articles/pi-chudnovsky/ 
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
#   If you want to add a Machin-like_formula, just add it to SET_OF_NAMES, SET_OF_MULTS etc.  Each list item Must contain the same
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
from datetime import timedelta
from functools import partial
import sys,time,multiprocessing,unicodedata,logging,os,argparse
try:
    # https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
    import Colorer
except ImportError:
    pass
if sys.version_info[0] < 3:
    print(os.path.basename(__file__) + " requires at least Python 3")
    sys.exit(1)
try:
    from gmpy2 import mpz,isqrt,mpfr,atan2,sqrt,get_context,const_pi  # Gumpy2 mpz large ints are ten times faster than python large int
except ImportError as exc:
    raise ImportError('This program requires gmpy2, please insatll. exiting....') from exc

# Change logging to INFO or WARNING to see less output
logging.basicConfig(level=("DEBUG"),format='[%(levelname)s] %(asctime)s %(funcName)s: %(processName)s %(message)s')

if sys.getdefaultencoding() != 'utf-8':
    logging.warning('char encoding is "{}" not utf-8 Some characters like pi may not print correctly'
            .format(sys.getdefaultencoding()) )

# CONSTANTS
LOG2_10 = 3.321928094887362
SAYPI = unicodedata.lookup("GREEK SMALL LETTER PI")
LAST_5_DIGITS_OF_PI = {
             10 : "26535",
            100 : "70679",
           1000 : "01989",
          10000 : "75678",
         100000 : "24646",
        1000000 : "58151",
        1234567 : "14707",
       10000000 : "55897",
      100000000 : "51592",
     1000000000 : "45519",
    }
#  Took values from lists from Machin and Miachin like formulae here:
#  https://en.wikipedia.org/wiki/Machin-like_formula
#  Add to all 4 lists to add a new formula  Each formula list entry *Must* be the same size for each type of list
SET_OF_NAMES = ["John Machin 1706",
    "F. C. M. Störmer 1896",
    "Kikuo Takano 1982",
    "Hwang Chien-Lih, 1997",
    "Hwang Chien-Lih, 2003",
    "Jörg Uwe Arndt 1993 ",
    "Hwang Chien-Lih, 2004",
    "\tRadius Generator- Fabrice Bellard?, 1997 \n\tπ = 126N∑n=0(−1)n210n(−254n+1−14n+3+2810n+1−2610n+3−2210n+5−2210n+7+110n+9)\n",
    "\tThe Square AGM - Salamin & Brent, 1976\n\tπ = limit as n goes to infinity  (an+bn)**2/(4tn)\n",
     "\tChudnovsky brothers  1988 \n\tπ = (Q(0, N) / 12T(0, N) + 12AQ(0, N))**(C**(3/2))\n",
     "\tconst_pi() function from the gmpy2 library" ]
SET_OF_MULTS = [ [4,1],  # Machin 1706
    [44,7,12,24], #  F. C. M. Störmer 1896
    [12,32,5,12], # Kikuo Takano (1982)
    [183,32,68,12,12,100], # Hwang 1997
    [183,32,68,12,100,12,12], # Hwang 2003
    [36462,135908,274509,39581,178477,114569,146571,61914,69044,89431,43938], # Jörg Uwe Arndt 1993
    [36462,26522,19275,3119,3833,5183,37185,11010,3880,16507,7476],  # Hwang 2004
    ["Place","holder"], # Improved Hex - Fabrice Bellard, 1997
    ["Place","holder"], # The Square AGM - Salamin & Brent, 1976
    ["Place","holder"],  # Chudnovsky brothers  1988
    ["Place","holder"] ]  # gmpy2 const_pi() 
SET_OF_DENOMS = [ [5,239],
    [57,239,682,12943],
    [49,57,239,110443],
    [239,1023,5832,110443,4841182,6826318],
    [239,1023,5832,113021,6826318,33366019650,43599522992503626068],
    [390112,485298,683982,1984933,2478328,3449051,18975991,22709274,24208144,201229582,2189376182],
    [51387,485298,683982,1984933,2478328,3449051,18975991,22709274,24208144,201229582,2189376182],
    ["Place","holder"], # Improved Hex - Fabrice Bellard, 1997
    ["Place","holder"], # The Square AGM - Salamin & Brent, 1976
    ["Place","holder"], # Chudnovsky brothers  1988
    ["Place","holder"] ]  # gmpy2 const_pi() 
SET_OF_OPERS = [ [1,-1],
    [1,1,-1,1],
    [1,1,-1,1],
    [1,1,-1,1,-1,-1],
    [1,1,-1,1,-1,-1,1],
    [1,1,1,-1,1,-1,-1,1,-1,-1,-1],
    [1,1,1,-1,-1,-1,-1,-1,1,-1,-1],
    ["Place","holder"], # Improved Hex/Radius - Fabrice Bellard, 1997
    ["Place","holder"], # The Square AGM - Salamin & Brent, 1976
    ["Place","holder"], # Chudnovsky brothers  1988
    ["Place","holder"] ]  # gmpy2 const_pi() 

NUM_OF_FORMULAE = len(SET_OF_DENOMS)
FROM_RANGE = "[1 to {}]".format(NUM_OF_FORMULAE)

# utility functions
def say_formula(credit,mults,denoms,signs):
    """ say formula creates a printed version of a formula using parts of the formula
        :param string credit: nmae and date of formula author
        :param list mults: array of Manchin multiples
        :param list denoms: array of arctan denominators
        :param list signs: array of plus or minus for formula arctans
        :return string: printed version of the formula
        """
    # format a Manchin like formula string from 3 lists and an author's credit or just print others
    if ('Chudnovsky' in credit) or ('AGM' in credit) or ("Bellard" in credit) or ("const_pi" in credit):
        form = credit
    else:
        form = "\t{}\n\t ".format(credit) + SAYPI + "/4 = "
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

def range_type(test_value, rngMin=1, rngMax=10):
    """ range_type is a simple range check for low high boundry
        :param string value: string to convert to a positive int (drop commas and other chars)
        :param int rngMin:  integer minimum allowed value
        :param int rngMax:  interger max allowed value
        :return int accepted value
        :throws argparse.ArgumentTypeError
        """
    value = int(''.join(filter(lambda i: i.isdigit(), test_value)) )  # chuck everthing but digits.  result is a positive int
    if value not in range(rngMin,rngMax+1):
        raise argparse.ArgumentTypeError('value {} not in range {:,} to {:,}'
            .format(value,rngMin,rngMax))
    return value

#Classes for various Pi formulae

class PiAGM:
    def __init__(self,ndigits):
        self.ndigits = ndigits
        self.cdigits = self.ndigits + len(str(self.ndigits))+9      # Extra digits to reduce trailing error More factors means more error
        self.iters = 0
        self.start_time = 0

    def compute(self):
        # Found formula here: https://www.kurims.kyoto-u.ac.jp/~ooura/pi_fft.html
        # This is an FFT modified AGM routine  POW() is not used 
        get_context().precision=int(self.cdigits * LOG2_10)
        epsilon = mpfr(1)/pow(mpfr(10),self.ndigits)
        logging.debug('AGM precision ({:,}) Started '
            .format(self.ndigits ) )
        self.start_time = time.time()
        c = mpfr(sqrt(0.125))
        a  = mpfr(1 + 3 * c )
        b  = mpfr(  sqrt(a))
        e = mpfr(b - 0.625)
        b *= 2  
        c = e - c
        a +=  e
        npow = 4
        
        while e > epsilon:
            npow *= 2
            e = (a + b) / 2
            b = sqrt(a * b)
            e = e - b
            b *= 2
            c = c - e
            a = b + e
            self.iters += 1
            if self.iters % 10  == 0:
                logging.debug('AGM ... {:,} iterations and {:.2f} seconds.'
                    .format( int(self.iters),time.time() - self.start_time))
        # a and b have converged to the AGM
        e = e * e / 4
        a = a + b
        get_context().precision=int((self.ndigits+2 ) * LOG2_10)
        pi =  (a * a - e - e / 2) / (a * c - e) / npow
        logging.debug('AGM Done! {:,} iterations and {:.2f} seconds.'
            .format(self.iters,time.time() - self.start_time) )
        return  str(pi)[:-2],self.iters,time.time()-self.start_time

class PiBellard:
    def __init__(self,ndigits):
        self.ndigits = ndigits
        self.iters = 0
        self.start_time = 0
        self.iter_time = 0

    def compute(self):
        """
        http://en.wikipedia.org/wiki/Bellard%27s_formula
        https://en.wikipedia.org/wiki/Bailey%E2%80%93Borwein%E2%80%93Plouffe_formula
        """
        cdigits = self.ndigits+15
        get_context().precision=int(cdigits * LOG2_10) # Precision isn't digits  
        self.start_time = time.time()  # Start the clock for total time
        logging.debug('BBP precision({:,}) Started '
            .format(self.ndigits ) )
        if self.ndigits > 10000:        
            logging.warning("\nWARNING\nWARNING Will Robinson\nBellard is a generator and will take a very long time if digits is > 10k.\n")
        pi = mpfr(0)
        
        self.iter_time = time.time()
        for i in range( ndigits):
            k = mpfr(i)
            a = mpfr(1/(pow(16,k)))
            b = mpfr(4/(8*k+1))
            c = mpfr(2/(8*k+4))
            d = mpfr(1/(8*k+5) )
            e = mpfr(1/(8*k+6) )
            r = mpfr(a*(b-c-d-e))
            pi += r
            self.iters += 1
            if self.iters % 10000  == 0:
                logging.debug('Bellard ... {:,} iterations and {:.2f} seconds 10k iters took {:.2f}.'
                    .format( int(self.iters),time.time() - self.start_time, time.time() - self.iter_time) )
                self.iter_time = time.time()
            
        get_context().precision=int((self.ndigits+4) * LOG2_10) # Precision isn't digits  need some math
        pi = pi + 0
        logging.debug('Bellard Done! {:,} iterations and {:.2f} seconds.'
            .format(self.iters,time.time() - self.start_time) )
        return str(pi)[:self.ndigits+2],self.iters,time.time()-self.start_time


class PiMachin:

    def __init__(self,ndigits,name,denoms,mults,operators):
        """ Initialization
        :param int digits: digits of PI computation
        :param string name: name for credit on formula
        :param list  denoms a lis of ints for the denomiators
        :param list  mults: a list of ints as multipliyers for the Machin formula
        :param list operators: a list of 1 or -1 to cause addition or subtraction
        """
        self.name = name  # Not used.  Planned to use in logging
        self.ndigits = ndigits
        self.denoms = denoms
        self.mults = mults
        self.operators = operators
        self.xdigits = len(denoms)+7              # Extra digits to reduce trailing error More factors means more error
        self.start_time = 0

    def ArctanDenom(self,d):

        cdigits = self.ndigits+self.xdigits
        get_context().precision=int(cdigits * LOG2_10)
        # Calculates arctan(1/d) = 1/d - 1/(3*d^3) + 1/(5*d^5) - 1/(7*d^7) + ...
        logging.debug('arctan(1/%d) Started  ',d )
        total = mpfr(0)
        arc_start_time = time.time()  # Start the clock for this arctan calulation
        total = mpfr(atan2(mpfr(1),mpfr(d)))
        logging.debug('arctan(1/{}) Done!   {:.2f} seconds.'
            .format(int(d),time.time() - arc_start_time))
        return total,int(1) # I used to calulate arctan by hand.  Now I just use atan2() so just one iteration here
    #
    def compute(self):
        self.start_time = time.time()  # Start the clock for total time
        ndigits = self.ndigits
        cdigits = self.ndigits+self.xdigits
        logging.info("Starting Machin-Like formula to {:,} decimal places"
            .format(ndigits) )
        get_context().precision=int(cdigits * LOG2_10)
        logging.debug("Starting %d Pool threads to calculate arctan values.",
            len(self.denoms) )
        p =  multiprocessing.Pool(processes=(len(self.denoms))) # get some threads for our pool
        results=p.map(self.ArctanDenom, self.denoms) # one thread per arctan(1/xxxx)
        p.close()
        p.join()  # wait for them to finish
        # Now we have the arctan calculations from the pool threads in results[]
        # Apply chosen Formula to the results and calculate pi using mults and signs
        logging.debug ("Now multiplying and summing all arctan results")
        arctanSum = pi = mpfr(0)
        iters = 0
        for i, result in enumerate(results):
            iters += result[1]  # Keep track of this thread's iterations for later
            arctanSum += mpfr(mpfr(self.mults[i])*mpfr(result[0])*mpfr(self.operators[i])) # Add or subtract the product from the accumulated arctans
        get_context().precision= int((self.ndigits+2) * LOG2_10)
        pi = mpfr(4) * arctanSum # change pi/4 = x to pi = 4 * x
        # We calculated extra digits to compensate for roundoff error.
        # Chop off the extra digits now.
        return str(pi)[:self.ndigits+2],iters,time.time()-self.start_time

class slow_chudnovsky:
    A = mpz(13591409)
    B = mpz(545140134)
    C = mpz(640320)
    D = mpz(426880)
    E = mpz(10005)
    """Slower, easer to understand, version of Chudnovsky Bros without Binary Splitting 
        Not used here  Just left over from earlier testing
    """
    def __init__(self,ndigits):
        """ Initialization
        :param int ndigits: digits of PI computation
        """
        self.ndigits = ndigits
        # 20 extra digits fluff for lots of calculations
        self.scale = pow(mpz(10),mpz(ndigits+20))
        self.iters = mpz(0)
        self.M10K = mpz(10000)
        self.k = mpz(1)
        self.a_k = mpz(self.scale)
        self.a_sum = mpz(self.scale)
        self.b_sum = mpz(0)
        self.C_cubed_over_24 = pow(self.C,mpz(3)) // mpz(24)

    def compute(self):
        get_context().precision =  int(LOG2_10 * (self.ndigits + 20) )
        self.start_time = time.time()
        logging.debug("Starting {} formula to {:,} decimal places"
                .format(name,self.ndigits) )
        while True:
            self.iters += mpz(1)
            if self.iters % self.M10K  == mpz(0):
                logging.debug('slow_chudnovsky ... {:,} iterations and {:.2f} seconds.'
                    .format( int(self.iters),time.time() - self.start_time))
            self.a_k *= mpz(-(mpz(6)*self.k-mpz(5)) * (mpz(2)*self.k-mpz(1)) * (mpz(6)*self.k-mpz(1)))
            self.a_k //=  mpz((self.k*self.k*self.k) * self.C_cubed_over_24)
            self.a_sum += mpz(self.a_k)
            self.b_sum += mpz(self.k * self.a_k)
            self.k += mpz(1)
            if self.a_k == mpz(0):
                break
        # convert to a fraction 
        total_sum = mpfr(self.A * self.a_sum + self.B * self.b_sum)
        pi = ((self.D * sqrt(self.E)) / total_sum)*self.scale
        logging.debug('{} calulation Done! {:,} iterations and {:.2f} seconds.'
                .format( name, int(self.iters),time.time() - self.start_time))
        # Chop off extra 20 digits 
        return str(pi)[:-20],int(self.iters),time.time() - self.start_time


class PiChudnovsky:
    """Version of Chudnovsky Bros using Binary Splitting 
        So far this is the winner for fastest time to a million digits on my older intel i7
        https://gist.github.com/komasaru/c3f5227513e1692c8fba42fe337316bc started here.  Mine is about 15% faster
    """
    A = mpz(13591409)
    B = mpz(545140134)
    C = mpz(640320)
    D = mpz(426880)
    E = mpz(10005)
    C3_24  = pow(C, mpz(3)) // mpz(24)
    #DIGITS_PER_TERM = math.log(53360 ** 3) / math.log(10)  #=> 14.181647462725476
    DIGITS_PER_TERM = 14.181647462725476
    MMILL = mpz(1000000)

    def __init__(self,ndigits):
        """ Initialization
        :param int ndigits: digits of PI computation
        """
        self.ndigits = ndigits
        self.n      = mpz(self.ndigits // self.DIGITS_PER_TERM + 1)
        self.prec   = mpz((self.ndigits + 1) * LOG2_10)
        self.one_sq = pow(mpz(10),mpz(2 * ndigits))
        self.sqrt_c = isqrt(self.E * self.one_sq)
        self.iters  = mpz(0)
        self.start_time = 0

    def compute(self):
        """ Computation """
        try:
            self.start_time = time.time()
            logging.debug("Starting {} formula to {:,} decimal places"
                .format(name,ndigits) )
            __, q, t = self.__bs(mpz(0), self.n)  # p is just for recursion
            pi = (q * self.D * self.sqrt_c) // t
            logging.debug('{} calulation Done! {:,} iterations and {:.2f} seconds.'
                .format( name, int(self.iters),time.time() - self.start_time))
            get_context().precision= int((self.ndigits+10) * LOG2_10)
            pi_s = pi.digits() # gmpy's digits() returns a string of the mpz int
            pi_o = pi_s[:1] + "." + pi_s[1:]
            return pi_o,int(self.iters),time.time() - self.start_time
        except Exception as e:
            print (e.message, e.args)
            raise

    def __bs(self, a, b):
        """ PQT computation by BSA(= Binary Splitting Algorithm)
        :param int a: positive integer
        :param int b: positive integer
        :return list [int p_ab, int q_ab, int t_ab]
        """
        try:
            self.iters += mpz(1)
            if self.iters % self.MMILL  == mpz(0):
                logging.debug('Chudnovsky ... {:,} iterations and {:.2f} seconds.'
                    .format( int(self.iters),time.time() - self.start_time))
            if a + mpz(1) == b:
                if a == mpz(0):
                    p_ab = q_ab = mpz(1)
                else:
                    p_ab = mpz((mpz(6) * a - mpz(5)) * (mpz(2) * a - mpz(1)) * (mpz(6) * a - mpz(1)))
                    q_ab = pow(a,mpz(3)) * self.C3_24
                t_ab = p_ab * (self.A + self.B * a)
                if a & 1:
                    t_ab *= mpz(-1)
            else:
                m = (a + b) // mpz(2)
                p_am, q_am, t_am = self.__bs(a, m)
                p_mb, q_mb, t_mb = self.__bs(m, b)
                p_ab = p_am * p_mb
                q_ab = q_am * q_mb
                t_ab = q_mb * t_am + p_am * t_mb
            return [p_ab, q_ab, t_ab]
        except Exception as e:
            print (e.message, e.args)
            raise

class PiConstant:
    def __init__(self,ndigits):
        """ Initialization
        :param int ndigits: digits of PI computation
        """
        self.ndigits = ndigits
        self.iters = 1
        self.start_time = 0

    def compute(self):
        """ Computation """
        logging.debug("Starting {} formula to {:,} decimal places"
                .format(name,self.ndigits) )
        precn = int((self.ndigits + 2) * LOG2_10) 
        self.start_time = time.time()
        my_pi = str(const_pi(precn))[:-2]
        logging.debug('{} calulation Done! {:,} iterations and {:.2f} seconds.'
                .format( name, int(self.iters),time.time() - self.start_time))
        return my_pi,self.iters,time.time()-self.start_time 

# Main for running one of the classes and saving the output
if __name__ == '__main__':

    # create strings for each formula for help.  Add here and to the DESC_STRING if you added a formula to the lists above
    FORMULA_LIST = ""
    for i in range(NUM_OF_FORMULAE):
        FORMULA_LIST += " {} {} \n".format(i+1,say_formula(SET_OF_NAMES[i],SET_OF_MULTS[i],SET_OF_DENOMS[i],SET_OF_OPERS[i]))

    pgmName =  os.path.basename(sys.argv[0])
    DESC_STRING = """ {0} runs an algoritym from a list to calulate Pi to a number of decimal places
      Default: {0} --digits 100000 --file pi.txt --alog 4

      So -d 100,000,000 will take a while to finish, -d 1,000,000 very quickly
      A last 5 digit validity check is done on powers of ten (10,...100,000,000)
 eg.  {0} --file elbow.txt -d 1000000 -a 10
      {0} -f test.txt -d 123,456

      List of Formulae:

{1} """.format(pgmName,FORMULA_LIST)
    parser = argparse.ArgumentParser(description=DESC_STRING,formatter_class=argparse.RawDescriptionHelpFormatter )
    # add expected arguments
    parser.add_argument('-f','--file', nargs='?', dest='filename', default='pi.txt',
                required=False,  help="File Name to write Pi to.. Default is [%(default)s]")
    parser.add_argument('-d','--digits', nargs=1, dest='max_digits', metavar="[1 to 1,000,000,000]", default=[100000],
                type=partial(range_type, rngMin=1, rngMax=1000000000), required=False,
                help="How many digits to calculate.  Default is %(default)s ")
    parser.add_argument('-a','--algo',nargs=1, dest='algo', metavar=FROM_RANGE, default=[4],
                type=partial(range_type, rngMin=1, rngMax=NUM_OF_FORMULAE), required=False, help="Which Machin(like) formula. Default is %(default)s")

    args = parser.parse_args(sys.argv[1:])
    if args.max_digits:
        ndigits = int(args.max_digits[0])
    if args.filename:
        outFileName =  args.filename
    if args.algo:
        algox =  args.algo[0] - 1

    start_time = time.time()  # Start the clock for total time

    name = SET_OF_NAMES[algox]  # pull the chosen formula list from the list of formulae
    denoms = SET_OF_DENOMS[algox]
    mults = SET_OF_MULTS[algox]
    operators = SET_OF_OPERS[algox]

    logging.info("Computing {} to {:,} digits."
            .format(SAYPI,ndigits))

    if 'Chudnovsky' in name:
        obj = PiChudnovsky(ndigits)
    else:
        if 'AGM' in name:
            obj = PiAGM(ndigits)
        else:
            if 'Bellard' in name:
                obj = PiBellard(ndigits)
            else:
                if "const_pi" in name:
                    obj = PiConstant(ndigits)
                else:
                    obj = PiMachin(ndigits,name,denoms,mults,operators)
    # Calculate Pi using selected formula
    pi,iters,time_to_calc = obj.compute()

    if ndigits in LAST_5_DIGITS_OF_PI:
        endDigits = pi[-5:]  # Pull the last 5 digits for a cross check
        if LAST_5_DIGITS_OF_PI[ndigits] == endDigits:
            logging.info("Last 5 digits of {} were {} as expected at offset {:,}"
                .format(SAYPI,endDigits,ndigits-5 ))
        else:
            logging.warning("\n\nWRONG WRONG WRONG\nLast 5 digits were {} and are WRONG should be {}\nWRONG WRONG WRONG\n"
                .format(endDigits,LAST_5_DIGITS_OF_PI[ndigits]) )
    else:
        logging.info("Did not check last 5 digits of result, {:,} wasn't in the list of known values"
            .format(ndigits) )

    startWrite = time.time()
    with open(outFileName, mode='wt',encoding="utf-8") as outfile:
        outfile.write(pi)
    time_to_write = time.time() - startWrite
    logging.info("Calculated {} to {:,} digits using a formula of:\n {} {} "
        .format(SAYPI,ndigits,algox+1,say_formula(name,mults,denoms,operators) ) )
    logging.debug('Wrote {:,} digits of {} to file {} in {}'
        .format(ndigits,SAYPI,outFileName,str(timedelta(seconds=time_to_write))))
    logging.info("Calculation took {:,} iterations and {}."
        .format(int(iters),str(timedelta(seconds=time_to_calc))) )
    sys.exit(0)
