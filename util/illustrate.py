
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot( (1970,2009), (0,1) , 'b', label='Observed')
ax.plot( (2009,2050), (1,2), 'b--', label='Current Trend')
ax.plot( (2009,2050), (1,3) , 'r', label='Future Scenario')

ax.set_title("Fictitious Climate Time Series")
ax.set_ylabel("Some variable")
ax.set_xlabel("Year") 

ax.legend(loc=9)

fig.savefig('delta_method.png')
