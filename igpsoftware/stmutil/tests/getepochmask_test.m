

epochs = 2010:.5:2021

epochmask = getepochmask(epochs,[ ])
epochmask = getepochmask(epochs,[ -Inf +Inf])
epochmask = getepochmask(epochs,[ -Inf 2012])
epochmask = getepochmask(epochs,[2012 2014])
epochmask = getepochmask(epochs,[2011 2013 ; 2018 2020 ])
epochmask = getepochmask(epochs,[2011 2013 ; 2018 2020 ; 2021 +Inf])

