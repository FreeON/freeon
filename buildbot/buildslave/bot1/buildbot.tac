
from twisted.application import service
from buildbot.slave.bot import BuildSlave

basedir = r'/home/buildslave/bot1'
buildmaster_host = 'rust'
port = 9989
slavename = 'bot1'
passwd = 'bot1passwd'
keepalive = 600
usepty = 1
umask = None

application = service.Application('buildslave')
s = BuildSlave(buildmaster_host, port, slavename, passwd, basedir,
               keepalive, usepty, umask=umask)
s.setServiceParent(application)

