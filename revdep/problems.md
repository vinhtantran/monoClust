# puls

<details>

* Version: 0.1.1
* GitHub: https://github.com/vinhtantran/puls
* Source code: https://github.com/cran/puls
* Date/Publication: 2020-12-09 08:20:03 UTC
* Number of recursive dependencies: 147

Run `revdep_details(, "puls")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
     ERROR
    Running the tests in 'tests/testthat.R' failed.
    Last 13 lines of output:
      
      == Warnings ====================================================================
      -- Warning (test-as_monoclust.puls.R:2:3): print error when coercing wrong PULS class to MonoClust class --
      as_MonoClust does not know how to handle object of class list.
      Backtrace:
       1. testthat::expect_warning(...) test-as_monoclust.puls.R:2:2
       8. monoClust:::as_MonoClust.default(a)
      
      == Failed tests ================================================================
      -- Failure (test-as_monoclust.puls.R:2:3): print error when coercing wrong PULS class to MonoClust class --
      `{ ... }` did not throw the expected warning.
      
      [ FAIL 1 | WARN 1 | SKIP 6 | PASS 16 ]
      Error: Test failures
      Execution halted
    ```

