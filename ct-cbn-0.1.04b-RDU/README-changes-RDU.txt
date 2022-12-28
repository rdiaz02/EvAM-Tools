This includes a bug fix by RDU for an unjustified stopping when (loglik_new < loglik). See my exchange of emails with Niko Beerenwinkel on 2017-06-09 and following days. Niko alerted, though, to a potential nonidentifiability issue. 

This version also outputs the initial lambda and lambdas from each iteration (which allows to get the lambda from the final iteration) as well as the likelihood at each iteration.

I keep the original file as ct-cbn.h-original


In the past, I have used, from R, the lambda_i, for a lambda from final iteration. But I finally used the lambda for the re-run. (The R code, which uses this code, though, expects a "lambda_i" to exist. But that can be removed without harm (need to modify some extra code, though, where I do things like "final_lambda") 

