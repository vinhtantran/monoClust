text.MonoClust <- function(x,abbrev,which,rel.loc.x,...){
    ## Set up some defaults and abbreviate, then use
    ## text.rpart.
	if(missing(which))which <-3
    # REMOVE: Tan, 3/1/15, remove intertia line
	# text(.62,.5,"Inertia Explained", srt=90)
	text.rpart(x,which=which,abbrev=abbrev,rel.loc.x=rel.loc.x,...)

	if(abbrev=='L'){
			vars  <- x$frame$var
			uvars <-unique(vars)
			names <- uvars[ uvars != "<leaf>"]
			nums <- paste("V",1:length(names),sep="")
			legend(mean(max(x$frame$loc),min(x$frame$loc)),.9,paste(nums,names,sep=" = "),bty='n')
	}
}
# test
