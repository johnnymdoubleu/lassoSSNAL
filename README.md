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
        <li><a href="#installation">Code Sample</a></li>
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

Lasso (Least Absolute Shrinkage and Selection Operator) is a popular and widely used regression method that does linear regression and variable selection in the field of statistics and machine learning. We initiated this project as final year project at University of Edinbrugh. In this project, we aim to present an inexact Semismooth Newton augmented Lagrangian method to solve Lasso problems. With the advantage of the second order sparsity of the problem, we can reduce the expensive computational cost and yield a fast and highly efficient algorithm that was first developed by [(Li et al. 2018)](https://arxiv.org/abs/1607.05428).

The algorithm was first developed in **Matlab** and named [Suitelasso](https://github.com/MatOpt/SuiteLasso). We benchmarked the original code and built an R package to aid the majority of statisticians who actively use **R** in their research.

Currently, the project has not been published to CRAN at the moment but we aim to publish in the near future.

<p align="right">(<a href="#top">back to top</a>)</p>



### Built With

This section should list any major frameworks/libraries used to bootstrap your project. Leave any add-ons/plugins for the acknowledgements section. Here are a few examples.
* [R](https://www.r-project.org/)
* [C++](https://docs.microsoft.com/en-us/cpp/cpp/?view=msvc-170)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites

This is an example of how to list things you need to use the software and how to install them.
* npm
  ```sh
  npm install npm@latest -g
  ```

### Installation

_Below is an example of how you can instruct your audience on installing and setting up your app. This template doesn't rely on any external dependencies or services._

1. Get a free API Key at [https://example.com](https://example.com)
2. Clone the repo
   ```sh
   git clone https://github.com/your_username_/Project-Name.git
   ```
3. Install NPM packages
   ```sh
   npm install
   ```
4. Enter your API in `config.js`
   ```js
   const API_KEY = 'ENTER YOUR API';
   ```

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [x] Add Changelog
- [x] Publish to CRAN
- [ ] Include Elastic Net 

See the [open issues](https://github.com/johnnymdoubleu/lassoSSNAL/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Project Link: [https://github.com/johnnymdoubleu/lassoSSNAL](https://github.com/johnnymdoubleu/lassoSSNAL)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
The project was 


While constructing this repository and readme file. We have reviewed the sources below.
* [Choose an Open Source License](https://choosealicense.com)
* [Img Shields](https://shields.io)
* [Readme Template](https://github.com/othneildrew/Best-README-Template#top)

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
