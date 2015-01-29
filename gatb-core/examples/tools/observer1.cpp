//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We define some notification information class.
class MyEventInfo : public EventInfo
{
public:
    MyEventInfo (const std::string& message) : EventInfo(0), _message(message) {}
    const std::string& getMessage ()  { return _message; }
private:
    std::string _message;
};

// We define some Observer class.
class MyObserver : public IObserver
{
public:
    void update (EventInfo* evt, ISubject* subject)
    {
        MyEventInfo* info = dynamic_cast<MyEventInfo*> (evt);
        if (info != 0)  {  std::cout << "Receiving: " << info->getMessage() << std::endl;  }
    }
};

/********************************************************************************/
/*                Usage of the Observer/Subject class                           */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // we define a subject instance
    ISubject* subject = new Subject ();

    // we create a specific observer
    IObserver* observer = new MyObserver ();

    // we attach the observer to the subject
    subject->addObserver (observer);

    // the subject sends some notification => should be received by our observer
    subject->notify (new MyEventInfo ("Message that should be received"));

    // we detach the observer from the subject
    subject->removeObserver (observer);

    // the subject sends some notification => should not be received by our observer
    subject->notify (new MyEventInfo ("Message that should NOT be received"));
}
//! [snippet1]
