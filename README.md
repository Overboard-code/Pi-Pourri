# Pi-Pourri
Several formulae for calculating 100 million digits of Pi in 3 minutes, using python and GMPY2 

I wanted to see how long it would take to calulate pi to a million places to answer a kid's question  

I found a page https://medium.com/@cosinekitty/how-to-calculate-a-million-digits-of-pi-d62ce3db8f58  that had a program using Machin's formula from 1706:

<img src="https://render.githubusercontent.com/render/math?math={\frac {\pi }{4}}=4\arctan {\frac {1}{5}}-\arctan {\frac {1}{239}}">

 

The article also included a link to several Machin-like formulae:  https://en.wikipedia.org/wiki/Machin-like_formula
I added several to the code   I also used GMPY's mpfr() and mpz types to speed things along instead of just big integer support in python  The GMPY2 library has a wide range of high precision functions.  

I added 7 Machin formulae like:</b>

<img src="https://render.githubusercontent.com/render/math?math=%7B%5Cdisplaystyle%20%7B%5Cbegin%7Baligned%7D%7B%5Cfrac%20%7B%5Cpi%20%7D%7B4%7D%7D%3D%26%5C%3B183%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B239%7D%7D%2B32%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B1023%7D%7D-68%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B5832%7D%7D%5C%5C%26%2B12%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B110443%7D%7D-12%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B4841182%7D%7D-100%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B6826318%7D%7D%5C%5C%5Cend%7Baligned%7D%7D%7D%0A%20%20%20%20">

Then I added Chudnovsky:</b>

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%7D%0Aa%20%20%20%20%20%26%3D%20%5Csum%5E%5Cinfty_%7Bk%3D0%7D%20%5Cfrac%7B(-1)%5Ek%20(6k)!%7D%7B(3k)!(k!)%5E3%20640320%5E%7B3k%7D%7D%20%5C%5C%0A%20%20%20%20%20%20%26%3D%201%0A%20%20%20%20%20%20%20%20%20%20-%20%5Cfrac%7B6%5Ccdot5%5Ccdot4%7D%7B(1)%5E3%20640320%5E3%7D%0A%20%20%20%20%20%20%20%20%20%20%2B%20%5Cfrac%7B12%5Ccdot11%5Ccdot10%5Ccdot9%5Ccdot8%5Ccdot7%7D%7B(2%5Ccdot1)%5E3%20640320%5E6%7D%0A%20%20%20%20%20%20%20%20%20%20-%20%5Cfrac%7B18%5Ccdot17%5Ccdots13%7D%7B(3%5Ccdot2%5Ccdot1)%5E3%20640320%5E%7B9%7D%7D%0A%20%20%20%20%20%20%20%20%20%20%2B%20%5Ccdots%20%5C%5C%0Ab%20%20%20%20%20%26%3D%20%5Csum%5E%5Cinfty_%7Bk%3D0%7D%20%5Cfrac%7B(-1)%5Ek%20(6k)!k%7D%7B(3k)!(k!)%5E3%20640320%5E%7B3k%7D%7D%20%5C%5C%0A%5Cfrac%7B1%7D%7B%5Cpi%7D%20%26%3D%20%5Cfrac%7B13591409a%20%2B%20545140134b%7D%7B426880%20%5Csqrt%7B10005%7D%7D%20%5C%5C%0A%5Cpi%20%20%20%20%20%20%20%20%20%20%20%26%3D%20%5Cfrac%7B426880%20%5Csqrt%7B10005%7D%7D%7B13591409a%20%2B%20545140134b%7D%0A%5Cend%7Balign%7D">

Finaly I added Arithmatic Geometric Mean </b>

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/824a061756f72d84359eba13d2e8bfcda777f9f4">

Here is the help for the program type ```pipourri.py -h```  to see it
```
usage: pi-pourri.py [-h] [-f [FILENAME]] [-d [1 to 100,000,000]] [-a [1 to 9]]

 pi-pourri.py runs an algorthym from a list to calulate Pi to a number of decimal places
      Default: pi-pourri.py --digits 100000 --file pi.txt --alog 4

      So -d 100,000,000 will take a while to finish, -d 1,000,000 very quickly
      A last 5 digit check is done on powers of ten (10,...10000000,100000000)
 eg.  pi-pourri.py --file elbow.txt -d 1000000 -a 4
      pi-pourri.py -f test.txt -d 123,456

      List of Formulae:

 1      John Machin 1706
         π/4 =  4*arctan(1/5)
                - arctan(1/239)

 2      F. C. M. Störmer 1896
         π/4 =  44*arctan(1/57)
                + 7*arctan(1/239)
                - 12*arctan(1/682)
                + 24*arctan(1/12943)

 3      Kikuo Takano 1982
         π/4 =  12*arctan(1/49)
                + 32*arctan(1/57)
                - 5*arctan(1/239)
                + 12*arctan(1/110443)

 4      Hwang Chien-Lih, 1997
         π/4 =  183*arctan(1/239)
                + 32*arctan(1/1023)
                - 68*arctan(1/5832)
                + 12*arctan(1/110443)
                - 12*arctan(1/4841182)
                - 100*arctan(1/6826318)

 5      Hwang Chien-Lih, 2003
         π/4 =  183*arctan(1/239)
                + 32*arctan(1/1023)
                - 68*arctan(1/5832)
                + 12*arctan(1/113021)
                - 100*arctan(1/6826318)
                - 12*arctan(1/33366019650)
                + 12*arctan(1/43599522992503626068)

 6      Jörg Uwe Arndt 1993
         π/4 =  36462*arctan(1/390112)
                + 135908*arctan(1/485298)
                + 274509*arctan(1/683982)
                - 39581*arctan(1/1984933)
                + 178477*arctan(1/2478328)
                - 114569*arctan(1/3449051)
                - 146571*arctan(1/18975991)
                + 61914*arctan(1/22709274)
                - 69044*arctan(1/24208144)
                - 89431*arctan(1/201229582)
                - 43938*arctan(1/2189376182)

 7      Hwang Chien-Lih, 2004
         π/4 =  36462*arctan(1/51387)
                + 26522*arctan(1/485298)
                + 19275*arctan(1/683982)
                - 3119*arctan(1/1984933)
                - 3833*arctan(1/2478328)
                - 5183*arctan(1/3449051)
                - 37185*arctan(1/18975991)
                - 11010*arctan(1/22709274)
                + 3880*arctan(1/24208144)
                - 16507*arctan(1/201229582)
                - 7476*arctan(1/2189376182)

 8      The Square AGM - Salamin & Brent, 1976
        π = limit as n goes to infinity  (an+bn)**2/(4tn)

 9      Chudnovsky brothers  1988
        π = (Q(0, N) / 12T(0, N) + 12AQ(0, N))**(C**(3/2))



options:
  -h, --help            show this help message and exit
  -f [FILENAME], --file [FILENAME]
                        File Name to write Pi to.. Default is pi.txt
  -d [1 to 100,000,000], --digits [1 to 100,000,000]
                        How many digits to calculate. Default is [100000]
  -a [1 to 9], --algo [1 to 9]
                        Which Machin(like) formula. Default is [4]
```
