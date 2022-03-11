<div id="top"></div>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <!--
  <a href="https://github.com/johnnymdoubleu/lassoSSNAL">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a> -->

  <h3 align="center">SSNAL Implementation in R</h3>

  <p align="center">
    Semismooth Newton Augmented Lagrangian method implementation in R to solve Lasso problems
    <br />
    <a href="https://github.com/johnnymdoubleu/lassoSSNAL"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/johnnymdoubleu/lassoSSNAL/Best-README-Template">View Demo</a>
    ·
    <a href="https://github.com/johnnymdoubleu/lassoSSNAL/issues">Report Bug</a>
    ·
    <a href="https://github.com/johnnymdoubleu/lassoSSNAL/issues">Request Features</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites & Installation</a></li>
        <!--li><a href="#installation">Code Sample</a></li-->
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

Lasso (Least Absolute Shrinkage and Selection Operator) is a popular and widely used regression method that does linear regression and variable selection in the field of statistics and machine learning. As a final year student at University of Edinbrugh (Myung Won Lee and Michael Renfrew), we initiated this as our final year project supervised by Dr. Daniel Paulin. In this project, we aim to present an inexact Semismooth Newton augmented Lagrangian method to solve Lasso problems. With the advantage of the second order sparsity of the problem, we can reduce the expensive computational cost and yield a fast and highly efficient algorithm that was first developed by [(Li et al. 2018)](https://arxiv.org/abs/1607.05428).

The algorithm was first developed in **Matlab** and named [Suitelasso](https://github.com/MatOpt/SuiteLasso). We benchmarked the original software and built an R package to aid the majority of statisticians who actively use **R** in their research.

*Currently, the project has not been published to CRAN at the moment but we aim to publish in the near future.*

<p align="right">(<a href="#top">back to top</a>)</p>



### Built With

This package is mainly built in **R** with the aid of **C++*.
* [R](https://www.r-project.org/)
* [C++](https://docs.microsoft.com/en-us/cpp/cpp/?view=msvc-170)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites & Installation

This is an example of how to list things you need to use the software and how to install them.
* Required Packages to run the algorithm
  - [pracma](https://github.com/cran/pracma)
  - [rmatio](https://github.com/stewid/rmatio)
  - [Rcpp](https://github.com/RcppCore/Rcpp)
  - [RSpectra](https://github.com/yixuan/RSpectra)
  ```R
  install.packages("pracma")
  install.packages("rmatio")
  install.packages("Rcpp")
  install.packages("RSpectra")
  ```
<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

1. Clone the repo
   ```sh
   git clone https://github.com/johnnymdoubleu/lassoSSNAL.git
   ```
2. Run R GUI / RStudio in the directory lassoSSNAL
3. Install / Load required packages
   ```R
   library(rmatio) 
   library(Rcpp) 
   library(RSpectra)
   library(Matrix)
   ```
4. Source **R** and **C++** files
   ```R
   source("Classic_Lasso_SSNAL.R")
   source("Classic_Lasso_SSNAL_main.R")
   source("Classic_Lasso_SSNCG.R")
   source("proj_inf.R")
   source("linsyssolve.R")
   source("findstep.R")
   source("psqmry.R")
   source("matvec_ClassicLasso.R")
   source("findnnz.R")
   sourceCpp("mex_matrix_mult.cpp")
   sourceCpp("mexsigma_update_classic_Lasso_SSNAL.cpp")
   ```
5. Define matrix/vector containing response and explanatory vairable by loading data
   ```R
   data <- read.mat("UCIdata/abalone_scale_expanded7.mat")
   A <- data$A
   b <- data$b
   ```
6. Define necessary arguments for the function
   ```R
   eps <- 2.220446e-16 # Copy the MATLAB eps essentially
   n <- ncol(A)
   c <- 10^(-4) # 10^(-3)
   rho <- c * max(abs(t(t(b) %*% A)))
   
   # evaluate the Lipschitz condition
   lipfun <- function(b, A){return(t(t(A%*%b) %*% A))}
   eigs_AtA <- eigs_sym(lipfun, k = 1, n = n, args = A)
   
   # set running parameters for function
   opts <- c()
   opts$stoptol <- 1e-6
   opts$Lip <- eigs_AtA$values
   opts$Ascale <- 1
   opts$maxiter <- 1000
   ```
7. Execute the SSNAL algorithm with profiling
   ```R
   Rprof()
   clo <- Classic_Lasso_SSNAL(A, b, n, rho, opts)
   Rprof(NULL)
   summaryRprof()
   
   # portion of the output values
   cat("min(X) = ", clo$info$minx, "\n")
   cat("max(X) = ", clo$info$max, "\n")
   cat("nnz = ", findnnz(clo$info$x,0.999)$k, "\n")
   ```

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [x] Add Changelog
- [ ] Publish to CRAN
- [ ] Include Elastic Net 

See the [open issues](https://github.com/johnnymdoubleu/lassoSSNAL/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b r.ssnal/AmazingUpdate`)
3. Commit your Changes (`git commit -m 'Add some AmazingUpdate'`)
4. Push to the Branch (`git push origin feature/AmazingUpdata`)
5. Open a Pull Request

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3.0. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Project Link: [https://github.com/johnnymdoubleu/lassoSSNAL](https://github.com/johnnymdoubleu/lassoSSNAL)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
We would like to express token of appreciation to our superviosr Dr. Daniel Paulin and our loved ones.

While constructing this repository and readme file. We have reviewed the sources below.
* [Choose an Open Source License](https://choosealicense.com)
* [Img Shields](https://shields.io)
* [Readme Template](https://github.com/othneildrew/Best-README-Template#top)

```
@thesis{Lee, Renfrew, thsis, 2022:rssnal,
   author    = {Lee, MyungWon and Renfrew, Michael},
   title     = {Semismooth Newton Augmented Lagrangian Method with Implementation in R},
   month     = {March},
   year      = {2022},
   howpublished = {Undergraduate Honors Thesis},
   school    = {University of Edinburgh},
   address   = {Edinburgh, United Kingdon}
   }
```

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/johnnymdoubleu/lassoSSNAL?color=light
[contributors-url]: https://github.com/johnnymdoubleu/lassoSSNAL/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/johnnymdoubleu/lassoSSNAL
[forks-url]: https://github.com/johnnymdoubleu/lassoSSNAL/network/members
[stars-shield]: https://img.shields.io/github/stars/johnnymdoubleu/lassoSSNAL
[stars-url]: hhttps://github.com/johnnymdoubleu/lassoSSNAL/stargazers
[issues-shield]: https://img.shields.io/github/issues/johnnymdoubleu/lassoSSNAL
[issues-url]: https://github.com/johnnymdoubleu/lassoSSNAL/issues
[license-shield]: https://img.shields.io/github/license/johnnymdoubleu/lassoSSNAL
[license-url]: https://github.com/johnnymdoubleu/lassoSSNAL/blob/main/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/johnnymwlee/
[product-screenshot]: images/screenshot.png
