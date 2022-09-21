# Pi-Pourri
Several formulae for calculating 100 million digits of Pi in less than 5 minutes, using python and GMPY2 

I wanted to see how long it would take to calulate pi to a million places to answer a kid's question.  
Turns out a million is about .6 seconds.
You can use Chudnovsky (--alog 10) to calculate 1 billion digits in about 45 minutes. Others take around an hour.  Manchin formulas take longer, I only had patience for 100 million at 35 minutes, Manchin would take about 9 hours for a billion I think. Memory gets to be a big deal for Manchin 30 or forty GIG.  
```

I found a page https://medium.com/@cosinekitty/how-to-calculate-a-million-digits-of-pi-d62ce3db8f58  that had a program using Machin's formula from 1706:

<img src="https://render.githubusercontent.com/render/math?math={\frac {\pi }{4}}=4\arctan {\frac {1}{5}}-\arctan {\frac {1}{239}}">


The article also included a link to several Machin-like formulae:  https://en.wikipedia.org/wiki/Machin-like_formula
I added several to the code   I also used GMPY's mpfr() and mpz types to speed things along instead of just big integer support in python  The GMPY2 library has a wide range of high precision functions.  

I added 7 Machin formulae like:</b>

<img src="https://render.githubusercontent.com/render/math?math=%7B%5Cdisplaystyle%20%7B%5Cbegin%7Baligned%7D%7B%5Cfrac%20%7B%5Cpi%20%7D%7B4%7D%7D%3D%26%5C%3B183%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B239%7D%7D%2B32%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B1023%7D%7D-68%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B5832%7D%7D%5C%5C%26%2B12%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B110443%7D%7D-12%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B4841182%7D%7D-100%5Carctan%20%7B%5Cfrac%20%7B1%7D%7B6826318%7D%7D%5C%5C%5Cend%7Baligned%7D%7D%7D%0A%20%20%20%20">

Then I added Chudnovsky (The clear speed winner):</b>

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%7D%0Aa%20%20%20%20%20%26%3D%20%5Csum%5E%5Cinfty_%7Bk%3D0%7D%20%5Cfrac%7B(-1)%5Ek%20(6k)!%7D%7B(3k)!(k!)%5E3%20640320%5E%7B3k%7D%7D%20%5C%5C%0A%20%20%20%20%20%20%26%3D%201%0A%20%20%20%20%20%20%20%20%20%20-%20%5Cfrac%7B6%5Ccdot5%5Ccdot4%7D%7B(1)%5E3%20640320%5E3%7D%0A%20%20%20%20%20%20%20%20%20%20%2B%20%5Cfrac%7B12%5Ccdot11%5Ccdot10%5Ccdot9%5Ccdot8%5Ccdot7%7D%7B(2%5Ccdot1)%5E3%20640320%5E6%7D%0A%20%20%20%20%20%20%20%20%20%20-%20%5Cfrac%7B18%5Ccdot17%5Ccdots13%7D%7B(3%5Ccdot2%5Ccdot1)%5E3%20640320%5E%7B9%7D%7D%0A%20%20%20%20%20%20%20%20%20%20%2B%20%5Ccdots%20%5C%5C%0Ab%20%20%20%20%20%26%3D%20%5Csum%5E%5Cinfty_%7Bk%3D0%7D%20%5Cfrac%7B(-1)%5Ek%20(6k)!k%7D%7B(3k)!(k!)%5E3%20640320%5E%7B3k%7D%7D%20%5C%5C%0A%5Cfrac%7B1%7D%7B%5Cpi%7D%20%26%3D%20%5Cfrac%7B13591409a%20%2B%20545140134b%7D%7B426880%20%5Csqrt%7B10005%7D%7D%20%5C%5C%0A%5Cpi%20%20%20%20%20%20%20%20%20%20%20%26%3D%20%5Cfrac%7B426880%20%5Csqrt%7B10005%7D%7D%7B13591409a%20%2B%20545140134b%7D%0A%5Cend%7Balign%7D">

Finaly I added Arithmatic Geometric Mean </b>

<img src="https://wikimedia.org/api/rest_v1/media/math/render/svg/824a061756f72d84359eba13d2e8bfcda777f9f4">

Here is a sample output of about 3 minutes for 100 million digits:
```
python pi-pourri.py -d 100,000,000 -a 10
[INFO] 2022-09-25 15:59:14,361 <module>: MainProcess Computing π to 100,000,000 digits.
[DEBUG] 2022-09-25 15:59:24,276 compute: MainProcess Starting 	Chudnovsky brothers  1988 
	π = (Q(0, N) / 12T(0, N) + 12AQ(0, N))**(C**(3/2))
 formula to 100,000,000 decimal places
[DEBUG] 2022-09-25 15:59:29,742 __bs: MainProcess Chudnovsky ... 1,000,000 iterations and 5.47 seconds.
[DEBUG] 2022-09-25 15:59:36,543 __bs: MainProcess Chudnovsky ... 2,000,000 iterations and 12.27 seconds.
[DEBUG] 2022-09-25 15:59:42,042 __bs: MainProcess Chudnovsky ... 3,000,000 iterations and 17.77 seconds.
[DEBUG] 2022-09-25 15:59:51,825 __bs: MainProcess Chudnovsky ... 4,000,000 iterations and 27.55 seconds.
[DEBUG] 2022-09-25 15:59:57,436 __bs: MainProcess Chudnovsky ... 5,000,000 iterations and 33.16 seconds.
[DEBUG] 2022-09-25 16:00:04,441 __bs: MainProcess Chudnovsky ... 6,000,000 iterations and 40.17 seconds.
[DEBUG] 2022-09-25 16:00:10,074 __bs: MainProcess Chudnovsky ... 7,000,000 iterations and 45.80 seconds.
[DEBUG] 2022-09-25 16:00:26,668 __bs: MainProcess Chudnovsky ... 8,000,000 iterations and 62.39 seconds.
[DEBUG] 2022-09-25 16:00:33,767 __bs: MainProcess Chudnovsky ... 9,000,000 iterations and 69.49 seconds.
[DEBUG] 2022-09-25 16:00:39,648 __bs: MainProcess Chudnovsky ... 10,000,000 iterations and 75.37 seconds.
[DEBUG] 2022-09-25 16:00:49,520 __bs: MainProcess Chudnovsky ... 11,000,000 iterations and 85.24 seconds.
[DEBUG] 2022-09-25 16:00:55,810 __bs: MainProcess Chudnovsky ... 12,000,000 iterations and 91.53 seconds.
[DEBUG] 2022-09-25 16:01:02,876 __bs: MainProcess Chudnovsky ... 13,000,000 iterations and 98.60 seconds.
[DEBUG] 2022-09-25 16:01:08,859 __bs: MainProcess Chudnovsky ... 14,000,000 iterations and 104.58 seconds.
[DEBUG] 2022-09-25 16:01:51,402 compute: MainProcess 	Chudnovsky brothers  1988 
	π = (Q(0, N) / 12T(0, N) + 12AQ(0, N))**(C**(3/2))
 calulation Done! 14,102,733 iterations and 147.13 seconds.
[INFO] 2022-09-25 16:02:16,558 <module>: MainProcess Last 5 digits of π were 51592 as expected at offset 99,999,995
[INFO] 2022-09-25 16:02:17,075 <module>: MainProcess Calculated π to 100,000,000 digits using a formula of:
 10 	Chudnovsky brothers  1988 
	π = (Q(0, N) / 12T(0, N) + 12AQ(0, N))**(C**(3/2))
 
[DEBUG] 2022-09-25 16:02:17,075 <module>: MainProcess Wrote 100,000,000 digits of π to file pi.txt in 0:00:00.516872
[INFO] 2022-09-25 16:02:17,075 <module>: MainProcess Calculation took 14,102,733 iterations and 0:02:52.269288.

```


Here is the help for the program type ```python3 pi-pourri.py -h```  to see it
```
usage: pi-pourri.py [-h] [-f [FILENAME]] [-d [1 to 100,000,000]]
                    [-a [1 to 10]]

 pi-pourri.py runs an algoritym from a list to calulate Pi to a number of decimal places
      Default: pi-pourri.py --digits 100000 --file pi.txt --alog 4

      So -d 100,000,000 will take a while to finish, -d 1,000,000 very quickly
      A last 5 digit check is done on powers of ten (10,...100,000,000)
 eg.  pi-pourri.py --file elbow.txt -d 1000000 -a 10
      pi-pourri.py -f test.txt -d 123,456

      List of Formulae:

 1 	John Machin 1706
	 π/4 =  4*arctan(1/5)
 		- arctan(1/239)
 
 2 	F. C. M. Störmer 1896
	 π/4 =  44*arctan(1/57)
 		+ 7*arctan(1/239)
 		- 12*arctan(1/682)
 		+ 24*arctan(1/12943)
 
 3 	Kikuo Takano 1982
	 π/4 =  12*arctan(1/49)
 		+ 32*arctan(1/57)
 		- 5*arctan(1/239)
 		+ 12*arctan(1/110443)
 
 4 	Hwang Chien-Lih, 1997
	 π/4 =  183*arctan(1/239)
 		+ 32*arctan(1/1023)
 		- 68*arctan(1/5832)
 		+ 12*arctan(1/110443)
 		- 12*arctan(1/4841182)
 		- 100*arctan(1/6826318)
 
 5 	Hwang Chien-Lih, 2003
	 π/4 =  183*arctan(1/239)
 		+ 32*arctan(1/1023)
 		- 68*arctan(1/5832)
 		+ 12*arctan(1/113021)
 		- 100*arctan(1/6826318)
 		- 12*arctan(1/33366019650)
 		+ 12*arctan(1/43599522992503626068)
 
 6 	Jörg Uwe Arndt 1993 
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
 
 7 	Hwang Chien-Lih, 2004
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
 
 8 	Radius Generator- Fabrice Bellard?, 1997 
	π = 126N∑n=0(−1)n210n(−254n+1−14n+3+2810n+1−2610n+3−2210n+5−2210n+7+110n+9)
 
 9 	The Square AGM - Salamin & Brent, 1976
	π = limit as n goes to infinity  (an+bn)**2/(4tn)
 
 10 	Chudnovsky brothers  1988 
	π = (Q(0, N) / 12T(0, N) + 12AQ(0, N))**(C**(3/2))
 
 11 	const_pi() function from the gmpy2 library 
 

options:
  -h, --help            show this help message and exit
  -f [FILENAME], --file [FILENAME]
                        File Name to write Pi to.. Default is pi.txt
  -d [1 to 100,000,000], --digits [1 to 100,000,000]
                        How many digits to calculate. Default is [100000]
  -a [1 to 11], --algo [1 to 11]
                        Which Machin(like) formula. Default is [4]
```
