select count(distinct mid)
from targets t 
where score >= -0.7 and score<= -0.6