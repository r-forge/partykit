## transformation function returning input angle on scale [-pi, pi)
# to be used at the begining of circtree/circforest
# arguments 'start' and 'end' specify the scale on which the angle is given originally
angle_trans <- function(angle, start = NULL, end = NULL)
{
  # if arguments 'start' and 'end. are not handed over one of the following 4 cases is assumed:
  # (-pi, pi]
  # [0, 2*pi)
  # (-180, 180]
  # [0, 360)

  stopifnot(is.numeric(angle))
  angle <- na.omit(angle)
  stopifnot(length(angle) > 0)
  
  if(is.null(start) & is.null(end)) {
    
    if(min(angle) >= 0) {
      if(any(angle > 2 *pi)){ 
        response_range <- c(0,360)
        warning(sprintf("Range of response is guessed to be between %.2f and %.2f, specify 'response_range' otherwise.", 
          response_range[1], response_range[2]), call. = FALSE)
      } else {
        response_range <- c(0,2*pi)
        warning("Range of response is guessed to be between 0 and 2pi, specify 'response_range' otherwise.", 
          call. = FALSE)
      }
    } else {
      if(any(abs(angle) > pi)){
        response_range <- c(-180,180)
        warning(sprintf("Range of response is guessed to be between %.2f and %.2f, specify 'response_range' otherwise.", 
          response_range[1], response_range[2]), call. = FALSE)
      }  else {
        response_range <- c(-pi,pi)
        warning("Range of response is guessed to be between -pi and pi, specify 'response_range' otherwise.", 
          call. = FALSE)
      }
    }
 
  } else if(is.null(start) | is.null(end)) {

    stop("arguments 'start' and 'end' can only be used if both are defined")

  } else {

    stopifnot(is.numeric(start) & length(start) == 1)
    stopifnot(is.numeric(end) & length(end) == 1)
    stopifnot(start < end)
    
    # check input arguments 'start' and 'end'
    if(start < 0) {
      if(abs(start) != end) stop("for negative values the scale has to be symmetric around 0")
    } else {
      if(start != 0) stop("for non-negative values the scale has to start at 0")
    }
    response_range <- c(start, end)
    
  }

  stopifnot(all(angle >= response_range[1]) & all(angle <= response_range[2]))
  
  # transfer to a scale of length 2*pi
  angle <- angle/diff(response_range)*2*pi
  
  if(response_range[1] < 0) {
    angle[angle<0] <- angle[angle<0] + 2*pi
  }
  if(any(angle > 2*pi) | any(angle < 0 )) 
    warning(paste0("Response values are transformed to ", ifelse(any(is.null(c(start, end))), 
      "approximated ", "") , "circular range."))
  angle <- angle %% (2*pi)
  # values are on scale [0, 2*pi]
  # transform to [-pi,pi)
  angle[angle>pi] <- angle[angle>pi] - 2*pi
  
  attr(angle, "response_range") <- response_range

  return(angle)
}

## transformation function returning calculated angle on required scale [start, end)
# to be used at the end of circtree/circforest
# arguments 'start' and 'end' need to be handed over/received from angle_trans
# calculated angle is on scale [-pi, pi)
angle_retrans <- function(angle, start = NULL, end = NULL)
{
  
  if(any(is.null(c(start, end))) && !is.null(attr(angle, "response_range"))) {
    response_range <- attr(angle, "response_range")
  } else if(any(is.null(c(start, end)))) {
    stop("arguments 'start' and 'end' must be provided")
  } else {
    response_range <- c(start, end)
  }
  
  if(response_range[1] == 0) angle[angle<0] <- angle[angle<0] + 2*pi
  if(response_range[2] > 2*pi) angle <- angle/(2*pi) * diff(response_range)
  
  return(angle)
}


# test
if(FALSE){
  a <- runif(100, 0, 2*pi)
  at <- angle_trans(a)
  art <- angle_retrans(at$angle, start = at$start_angle, at$end_angle)
  summary(art)
  summary(a)
  any(a!=art)
  cbind(a, art, a-art, a== art)
  at$start_angle
  at$end_angle
  
  a <- seq(0, 24, by = 2)
  at <- angle_trans(a, 0, 24)
  art <- angle_retrans(at$angle, start = at$start_angle, at$end_angle)
  summary(art)
  summary(a)
  any(a!=art)
  cbind(a, art, a-art, a== art)
  at$start_angle
  at$end_angle
}
