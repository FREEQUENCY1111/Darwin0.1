"""Base Agent class and message protocol for the Darwin Council.

Every agent in Darwin follows this contract:
  - Has a unique name and role description
  - Communicates via AgentMessages through the shared context
  - Reports its status (idle / working / done / error)
  - Can vote on decisions when the council deliberates
  - Logs all actions for auditability

Architecture:
                    ┌──────────┐
         ┌─────────┤Processor ├─────────┐
         │         └──────────┘         │
    ┌────┴─────┐   ┌──────────┐   ┌────┴───────┐
    │Annotator ├───┤Orchestr. ├───┤Scrutinizer │
    └────┬─────┘   └──────────┘   └────┬───────┘
         │         ┌──────────┐         │
         └─────────┤Synthesiz.├─────────┘
                   └─────┬────┘
                    ┌────┴─────┐
                    │ Uploader │
                    └──────────┘

  Support layer: Context Jimi, Page Rational, Deep Research Bob, Code-checker Chris
"""

from __future__ import annotations

import time
import uuid
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, TYPE_CHECKING

from darwin.utils.logging import get_logger

if TYPE_CHECKING:
    from darwin.agents.context import ContextJimi


class MessageType(str, Enum):
    """Types of messages agents can exchange."""
    DATA = "data"              # passing data between agents
    REQUEST = "request"        # asking another agent to do something
    RESPONSE = "response"      # replying to a request
    VOTE = "vote"              # casting a vote in a council decision
    ALERT = "alert"            # flagging an issue
    STATUS = "status"          # status update
    LOG = "log"                # audit log entry


class AgentStatus(str, Enum):
    IDLE = "idle"
    WORKING = "working"
    DONE = "done"
    ERROR = "error"
    WAITING = "waiting"        # waiting on another agent


@dataclass
class AgentMessage:
    """A message passed between agents through the shared context."""
    sender: str
    recipient: str             # agent name or "council" for broadcast
    msg_type: MessageType
    payload: dict[str, Any]
    timestamp: float = field(default_factory=time.time)
    msg_id: str = field(default_factory=lambda: uuid.uuid4().hex[:8])

    def to_dict(self) -> dict:
        return {
            "msg_id": self.msg_id,
            "sender": self.sender,
            "recipient": self.recipient,
            "type": self.msg_type.value,
            "payload": self.payload,
            "timestamp": self.timestamp,
        }


class Agent(ABC):
    """Base class for all Darwin Council agents.

    Every agent:
      - Has a name and role
      - Connects to Context Jimi (shared state)
      - Can send/receive messages
      - Reports its status
      - Has an execute() method that does its work
    """

    name: str = "base_agent"
    role: str = "Base agent with no specialization"

    def __init__(self) -> None:
        self.log = get_logger(f"darwin.agent.{self.name}")
        self.status = AgentStatus.IDLE
        self._context: ContextJimi | None = None
        self._message_log: list[AgentMessage] = []

    def connect(self, context: ContextJimi) -> None:
        """Connect this agent to the shared context (Jimi)."""
        self._context = context
        context.register_agent(self)
        self.log.info(f"[bold]{self.name}[/] connected to context")

    @property
    def context(self) -> ContextJimi:
        if self._context is None:
            raise RuntimeError(f"{self.name} is not connected to Context Jimi")
        return self._context

    def send(self, recipient: str, msg_type: MessageType, payload: dict) -> AgentMessage:
        """Send a message to another agent through the context."""
        msg = AgentMessage(
            sender=self.name,
            recipient=recipient,
            msg_type=msg_type,
            payload=payload,
        )
        self._message_log.append(msg)
        self.context.route_message(msg)
        return msg

    def broadcast(self, msg_type: MessageType, payload: dict) -> AgentMessage:
        """Broadcast a message to the entire council."""
        return self.send("council", msg_type, payload)

    def receive(self, message: AgentMessage) -> None:
        """Called when a message is routed to this agent."""
        self._message_log.append(message)
        self.handle_message(message)

    def handle_message(self, message: AgentMessage) -> None:
        """Override to handle incoming messages. Default: log and ignore."""
        self.log.debug(
            f"  [{self.name}] received {message.msg_type.value} from {message.sender}"
        )

    def set_status(self, status: AgentStatus) -> None:
        """Update agent status and notify the context."""
        self.status = status
        self.context.update_agent_status(self.name, status)

    @abstractmethod
    def execute(self, **kwargs: Any) -> dict[str, Any]:
        """Execute this agent's primary task.

        Returns a result dict that will be stored in the shared context.
        """
        ...

    @abstractmethod
    def validate_input(self, **kwargs: Any) -> bool:
        """Check that required inputs are available before executing."""
        ...

    def get_audit_log(self) -> list[dict]:
        """Return this agent's complete message log."""
        return [m.to_dict() for m in self._message_log]
