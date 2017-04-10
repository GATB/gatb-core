#pragma once
// from http://stackoverflow.com/questions/36762248/why-is-stdqueue-not-thread-safe

#include <queue>
#include <mutex>
#include <condition_variable>

    template <typename T>
    class SharedVectorQueue
    {
    public:
        SharedVectorQueue(int nb_elts) : nb_elts(nb_elts) {
        
            std::vector<std::mutex>(/*nb_elts*/ 1).swap(mutex_);

            queue_.resize(nb_elts);
        };

        ~SharedVectorQueue(){
            for (int i = 0; i < nb_elts; i++)
            {
                //std::deque<T>().swap(queue_[i]); 
                std::vector<T>().swap(queue_[i]); 
            }
            std::vector<std::vector<T>>().swap(queue_); 
            std::vector<std::mutex>().swap(mutex_); 
        };

        SharedVectorQueue(const SharedVectorQueue&) = delete;

        bool pop_immediately(int i, T& item)
        {
            if (queue_[i].empty())
            {
                queue_[i].shrink_to_fit(); // prevents higher mem usage, huh
                return false;
            }
             // std::vector queue
            item = queue_[i].back();
            queue_[i].pop_back();
            
            /*
            item = queue_[i].front();
            queue_[i].pop_front();
            */
            return true;
        }

        void push_back(int i, const T& item)
        {
            std::unique_lock<std::mutex> mlock(mutex_[0]);
            queue_[i].push_back(item);
            mlock.unlock();     // unlock before notificiation to minimize mutex con

        }
        void push_back(int i, T&& item)
        {
            std::unique_lock<std::mutex> mlock(mutex_[0]);
            queue_[i].push_back(std::move(item));
            mlock.unlock();     // unlock before notificiation to minimize mutex con
        }


        int size(int i)
        {
            std::unique_lock<std::mutex> mlock(mutex_[0]);
            int size = queue_[i].size();
            mlock.unlock();
            return size;
        }    
        
        bool empty(int i)
        {
            return queue_[i].empty();
        }


    private:
        int nb_elts;
        std::vector<std::vector<T>> queue_;
        std::vector<std::mutex> mutex_;
        std::vector<std::condition_variable> cond_;
    }; 


         

    

    
