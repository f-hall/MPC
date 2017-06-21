I used the framework given by udacity for this project and took some code from the lessons for this one.

I implemeted the MPC to get the car to stay on track. The MPC is roughly just a standard implementation without much of special extras. I tested different constraints and furthermore tried quite a lot of variables to get a good result (that was the main work and I spend the most time for that). To avoid problems with the latency I used a clue from the forum and tried to implement it - it works really well - look into my code for this! Additionally I tested some speed caps and it seems to work fine until somewhere around 60 mph (where it begins to slinger a bit on the street, but still stays on track).

I commented the code, but if there are still questions just write me back in the review so I can correct the problems.

P.S: I used the lowest resolution and smallest window for the simulator, if that is important. There were some problems with the simulator for my T2P4, where everything was going well for me, but my reviewer had different results, where the car would leave the track. I think the simulator is still really early in development and there may be some difficulties to get same results on different Computer/OS - so I send a video, where you can see that it works for me as expected!

I hope everything works as expected for you, too.

Thanks for the review and have a nice week. :)

Greetings Frank


