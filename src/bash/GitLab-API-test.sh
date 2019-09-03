#! /bin/bash -l
# a query to retrieve all events
# could have been used to check all the users that joined, but the history is incomplete
# https://microasp.upsc.se/api/v4/projects/2/events?private_token=AH1gfzASjFjQTxyyhTqx&action=joined&per_page=200&after=2015-01-01

# retrieve all users, we go for 10 pages, a 100 entry each, just in case
for i in {1..10}; do

    # first we find the blocked users and delete them
    for uid in `curl -s "https://microasp.upsc.se/api/v4/users?private_token=ee7S_qa7_zgNCLjfJoMx&per_page=100&page=$i&blocked=true" | jq -r '.[] | {id} | tostring' | sed 's/.*://g' | sed 's:}::'`;
    do
	uname=`curl -s --header "PRIVATE-TOKEN: ee7S_qa7_zgNCLjfJoMx" "https://microasp.upsc.se/api/v4/users/$uid" | jq -r '.[] | {username}'`
	echo "Blocked user $uname ($uid) will be deleted"
	curl --request DELETE --header "PRIVATE-TOKEN: ee7S_qa7_zgNCLjfJoMx" "https://microasp.upsc.se/api/v4/users/$uid"
	if [ $? -ne 0 ]; then
	    echo "WARNING: User $uname ($uid) deletion failed"
	fi
    done

    # next we find the user without any recorded events and block them

    for uid in `curl -s "https://microasp.upsc.se/api/v4/users?private_token=ee7S_qa7_zgNCLjfJoMx&per_page=100&page=$i" | jq -r '.[] | {id} | tostring' | sed 's/.*://g' | sed 's:}::'`; do
    	events=`curl -s "https://microasp.upsc.se/api/v4/users/$uid/events?private_token=ee7S_qa7_zgNCLjfJoMx" | jq -r '.[]' | wc -l`
	if [ $events -eq 0 ]; then
	   uname=`curl -s "https://microasp.upsc.se/api/v4/users/$uid?private_token=ee7S_qa7_zgNCLjfJoMx" | jq -r '.[] | {username}'`
	   echo "User $uname ($uid) will be blocked"
	   curl --request POST --header "PRIVATE-TOKEN: ee7S_qa7_zgNCLjfJoMx" "https://microasp.upsc.se/api/v4/users/$uid/block"
	   if [ $? -ne 0 ]; then
	    echo "WARNING: User $uname ($uid) blocking failed"
	   fi
	fi 
    done
done



