## Release summary
A bug fix release with some changes in function argument and a warning message format.

## Test environments
* local Windows 10 installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.2
* macos-10.15.7 (on Github Actions), R 4.0.3
* win-builder (devel)
* R-hub windows-x86_64-devel (r-devel)
* R-hub ubuntu-gcc-release (r-release)
* R-hub fedora-clang-devel (r-devel)

## R CMD check results

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded

## Further comments

Fixed the lifecycle URL in README.md to avoid redirect.

## Reverse Dependency

There is one package, puls, that is reverse dependency on this package. It will fail in one unit test because of a very small change (remove a space before the end point in the message) in the warning message in as_MonoClust() function. **It is totally expected and because I'm also the maintainer of puls, I will push a fix to that package very soon**.
